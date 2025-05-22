import torch
import torch.nn.functional as F
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch_geometric.nn import SAGEConv
from torch_geometric.utils import negative_sampling
from torch_geometric.data import Data
from sklearn.preprocessing import StandardScaler

RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/PPInetwork_Embedding_view/"


def load_data():
    """
    Load graph data: nodes, edges, and node features
    """
    node_dict = {}
    # Read protein nodes and assign each an index
    with open(RAWDATA_DIR + "PPInetwork_proteins.txt") as f:
        node_list = f.readlines()
        node_list = [node.split("\n")[0] for node in node_list]
        index = 0
        for n in node_list:
            node_dict[n] = index
            index += 1
    # Read edge list (protein interactions) and map node names to indices
    with open(RAWDATA_DIR + "PPInetwork_edgelist.csv") as f:
        edge_list = f.readlines()
        edge_list.pop(0)
        edge_list = [inter.split(",") for inter in edge_list]
        edge_list = [(u, v.split("\n")[0]) for u, v in edge_list]
        edge_list = [[node_dict[u], node_dict[v]] for u, v in edge_list]
        edge_list = torch.tensor(edge_list)

    # Load node features from topology file
    x = torch.zeros(size=(len(node_dict), 10))
    with open("../../data/PPInetwork_topology_view/" + "PPInetwork_topology.txt") as f:
        lines = f.readlines()
        for tab in lines[1:]:  # Skip header
            tab = tab.split("\t")
            tab[-1] = tab[-1].split("\n")[0]
            prot = tab[0]
            features = tab[1:]
            features = [float(i) for i in features]
            x[node_dict[prot]] = torch.tensor(features)

    # Normalize features using standard scaling
    scaler = StandardScaler()
    data = Data(x=x, edge_index=edge_list.t().contiguous())
    data.x = torch.tensor(scaler.fit_transform(data.x.numpy()), dtype=torch.float)

    return data, node_dict


class GraphSAGE(torch.nn.Module):
    """
    GraphSAGE model with two convolutional layers and dropout
    """
    def __init__(self, in_channels, hidden_channels, out_channels, dropout_rate=0.4):
        super().__init__()
        self.conv1 = SAGEConv(in_channels, hidden_channels)
        self.conv2 = SAGEConv(hidden_channels, out_channels)
        self.dropout = dropout_rate

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, p=self.dropout, training=self.training)
        x = self.conv2(x, edge_index)
        return x


# Contrastive loss based on InfoNCE to learn node embeddings
def contrastive_loss(embeddings, edge_index, temperature=0.2, num_negatives=2):
    pos_i, pos_j = edge_index
    z = F.normalize(embeddings, dim=1)

    # Positive pairs similarity (connected nodes)
    z_i = z[pos_i]
    z_j = z[pos_j]
    pos_sim = (z_i * z_j).sum(dim=-1) / temperature  # (B,)

    # Negative pairs similarity (randomly sampled edges)
    neg_edge_index = negative_sampling(
        edge_index,
        num_nodes=embeddings.size(0),
        num_neg_samples=pos_i.size(0) * num_negatives
    )
    neg_i, neg_j = neg_edge_index
    z_ni = z[neg_i]
    z_nj = z[neg_j]
    neg_sim = (z_ni * z_nj).sum(dim=-1) / temperature

    # Combine positives and negatives for InfoNCE loss
    all_sims = torch.cat([
        pos_sim.unsqueeze(1),
        neg_sim.view(pos_sim.size(0), num_negatives)
    ], dim=1)

    labels = torch.zeros(pos_sim.size(0), dtype=torch.long, device=embeddings.device)  # positives are label 0

    return F.cross_entropy(all_sims, labels)


def train(model, optimizer, x, edge_index, epochs=50):
    """
    Training loop for the model
    """
    scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.8, patience=5, threshold=1e-2)
    for epoch in range(epochs):
        optimizer.zero_grad()
        embeddings = model(x, edge_index)
        loss = contrastive_loss(embeddings, edge_index)
        loss.backward()
        optimizer.step()
        scheduler.step(loss.item())
        if epoch % 10 == 0:
            lr = optimizer.param_groups[0]['lr']
            print(f"Epoch {epoch}, Loss: {loss.item():.4f}, LR: {lr:.6f}")
            # print(optimizer.state_dict)


def save_embs(model, data, device='cpu'):
    """
    Save learned node embeddings to a text file
    """
    model.eval()
    model.to(device)
    data = data.to(device)
    with torch.no_grad():
        embeddings = model(data.x, data.edge_index).cpu()
    node_dict = {}
    with open(RAWDATA_DIR + "PPInetwork_proteins.txt") as f:
        node_list = f.readlines()
        node_list = [node.split("\n")[0] for node in node_list]
        for index, n in enumerate(node_list):
            node_dict[index] = n
    # Write embeddings line-by-line with node name followed by embedding vector
    with open(VIEW_DIR + "SAGE_PPInetwork.txt", "w") as f:
        for i in range(embeddings.size(0)):
            node_name = node_dict[i]  # Récupère le nom du nœud
            embedding_values = " ".join(map(str, embeddings[i].tolist()))  # Convertit les embeddings en string
            f.write(f"{node_name} {embedding_values}\n")


if __name__ == '__main__':
    # Change device to 'cuda' if GPU is available
    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    device = "cpu"
    data, node_dict = load_data()
    model = GraphSAGE(in_channels=data.num_node_features, hidden_channels=32, out_channels=128)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    train(model, optimizer, data.x, data.edge_index, epochs=300)
    save_embs(model, data)
