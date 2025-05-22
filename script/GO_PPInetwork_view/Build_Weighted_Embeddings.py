import numpy as np
import torch
import torch.nn.functional as F
from torch.nn import Dropout, ReLU, Linear
from torch_geometric.nn import SAGEConv, MessagePassing
from torch_geometric.utils import negative_sampling
from torch_geometric.data import Data
from torch.optim.lr_scheduler import ReduceLROnPlateau
from sklearn.preprocessing import StandardScaler

RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/GO_PPInetwork_view/"


def load_data(edge_attr_file):
    """
    Load graph data: nodes, edges, and edge features
    """
    node_dict = {}
    # Read protein nodes and assign each an index
    with open(RAWDATA_DIR + "PPInetwork_proteins.txt") as f:
        node_list = [line.strip() for line in f]
        for index, n in enumerate(node_list):
            node_dict[n] = index
    # Read edge list (protein interactions) and map node names to indices
    with open(RAWDATA_DIR + "PPInetwork_edgelist.csv") as f:
        next(f)
        edge_index = []
        for line in f:
            u, v = line.strip().split(",")
            edge_index.append([node_dict[u], node_dict[v]])
    edge_index = torch.tensor(edge_index).t().contiguous()

    # Read edge attributes from given file
    edge_attr_dict = {}
    with open(VIEW_DIR + edge_attr_file) as f:
        for line in f:
            u, v, val = line.strip().split()
            uid, vid = node_dict[u], node_dict[v]
            edge_attr_dict[(uid, vid)] = float(val)
            edge_attr_dict[(vid, uid)] = float(val)

    # Build edge attribute tensor aligned with edge_index order
    edge_attr = []
    for i in range(edge_index.size(1)):
        src, tgt = edge_index[0, i].item(), edge_index[1, i].item()
        edge_attr.append(edge_attr_dict.get((src, tgt), 0.0))
    # print(np.mean(edge_attr), np.std(edge_attr))
    edge_attr = torch.tensor(edge_attr, dtype=torch.float).unsqueeze(1)

    # Fake node features (all zeros, dimension 1)
    x = torch.zeros(size=(len(node_dict), 1))

    scaler = StandardScaler()
    edge_attr = scaler.fit_transform(edge_attr.numpy())
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    return data, node_dict


class EdgeAwareSAGEConv(MessagePassing):
    """
    Custom GraphSAGE convolution layer that incorporates edge attributes
    """
    def __init__(self, in_channels, edge_dim, out_channels):
        super().__init__(aggr='mean')
        self.lin_node = Linear(in_channels, out_channels)
        self.lin_edge = Linear(edge_dim, out_channels)
        self.lin_update = Linear(in_channels + out_channels, out_channels)

    def forward(self, x, edge_index, edge_attr):
        return self.propagate(edge_index, x=x, edge_attr=edge_attr)

    def message(self, x_j, edge_attr):
        return self.lin_node(x_j) + self.lin_edge(edge_attr)

    def update(self, aggr_out, x):
        h = torch.cat([x, aggr_out], dim=1)
        return F.relu(self.lin_update(h))


class GraphSAGE(torch.nn.Module):
    """
    GraphSAGE model using two edge-aware convolution layers, with layer norm and dropout
    """
    def __init__(self, in_channels, edge_dim, hidden_channels, out_channels, dropout=0.4):
        super().__init__()
        self.conv1 = EdgeAwareSAGEConv(in_channels, edge_dim, hidden_channels)
        self.conv2 = EdgeAwareSAGEConv(hidden_channels, edge_dim, out_channels)
        self.norm1 = torch.nn.LayerNorm(hidden_channels)
        self.norm2 = torch.nn.LayerNorm(out_channels)

        self.dropout = dropout

    def forward(self, x, edge_index, edge_attr):
        x = self.norm1(self.conv1(x, edge_index, edge_attr))
        x = F.dropout(x, p=self.dropout, training=self.training)
        x = self.norm2(self.conv2(x, edge_index, edge_attr))
        return x


def contrastive_loss(embeddings, edge_index, temperature=0.2, num_negatives=2):
    """
    Contrastive InfoNCE loss for training embeddings using positive and negative edges
    """
    pos_i, pos_j = edge_index
    z = F.normalize(embeddings, dim=1)

    # Positive pairs
    z_i = z[pos_i]
    z_j = z[pos_j]
    pos_sim = (z_i * z_j).sum(dim=-1) / temperature  # (B,)

    # Negative pairs
    neg_edge_index = negative_sampling(
        edge_index,
        num_nodes=embeddings.size(0),
        num_neg_samples=pos_i.size(0) * num_negatives
    )
    neg_i, neg_j = neg_edge_index
    z_ni = z[neg_i]
    z_nj = z[neg_j]
    neg_sim = (z_ni * z_nj).sum(dim=-1) / temperature  # (B * num_negatives,)

    # Combine positives and negatives for InfoNCE loss
    all_sims = torch.cat([
        pos_sim.unsqueeze(1),  # (B, 1)
        neg_sim.view(pos_sim.size(0), num_negatives)  # (B, N)
    ], dim=1)  # (B, N+1)

    labels = torch.zeros(pos_sim.size(0), dtype=torch.long, device=embeddings.device)  # positives à l'indice 0

    return F.cross_entropy(all_sims, labels)


def train(model, optimizer, data, epochs=100, device="cpu"):
    """
    Training loop for the model
    """
    model.to(device)
    data = data.to(device)

    scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.9, patience=10, threshold=1e-2)

    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        embeddings = model(data.x, data.edge_index, data.edge_attr)
        loss = contrastive_loss(embeddings, data.edge_index)
        loss.backward()
        optimizer.step()

        scheduler.step(loss.item())

        if epoch % 10 == 0:
            lr = optimizer.param_groups[0]['lr']
            print(f"Epoch {epoch}, Loss: {loss.item():.4f}, LR: {lr:.6f}")


def save_embs(model, data, out_file):
    """
    Save learned node embeddings to a text file
    """
    model.eval()
    with torch.no_grad():
        embeddings = model(data.x, data.edge_index, data.edge_attr)
    with open(RAWDATA_DIR + "PPInetwork_proteins.txt") as f:
        node_list = [line.strip() for line in f]
    with open(VIEW_DIR + out_file, "w") as f:
        for i in range(embeddings.size(0)):
            f.write(f"{node_list[i]} {' '.join(map(str, embeddings[i].tolist()))}\n")


def run_pipeline(filename, embed_output_name, hidden_dim=16, out_dim=128, epochs=500, lr=0.01):
    """
    Main pipeline to load data, train model, and save embeddings
    """
    data, _ = load_data(filename)
    model = GraphSAGE(in_channels=data.num_node_features, edge_dim=1,
                      hidden_channels=hidden_dim, out_channels=out_dim)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    train(model, optimizer, data, epochs=epochs)
    save_embs(model, data, embed_output_name)


if __name__ == '__main__':
    # Run two pipelines with different GO term-based edge attribute files
    run_pipeline("GO-BP_PPInetwork.txt", "GO-BP_PPInetwork_embed.txt")
    run_pipeline("GO-CC_PPInetwork.txt", "GO-CC_PPInetwork_embed.txt")
