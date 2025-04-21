import torch
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv
from torch_geometric.utils import negative_sampling
from torch_geometric.data import Data
from sklearn.preprocessing import StandardScaler

RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/GO_PPInetwork_view/"


def load_data(edge_attr_file):
    node_dict = {}
    with open(RAWDATA_DIR + "PPInetwork_proteins.txt") as f:
        node_list = [line.strip() for line in f]
        for index, n in enumerate(node_list):
            node_dict[n] = index

    # Chargement des arêtes
    with open(RAWDATA_DIR + "PPInetwork_edgelist.csv") as f:
        edge_list = f.readlines()[1:]  # skip header
        edge_list = [line.strip().split(",") for line in edge_list]
        edge_index = [[node_dict[u], node_dict[v]] for u, v in edge_list]
        edge_index = torch.tensor(edge_index).t().contiguous()

    # Chargement des features de nœuds
    x = torch.zeros(size=(len(node_dict), 1))
    """with open(VIEW_DIR + "PPInetwork_topology.txt") as f:
        for line in f.readlines()[1:]:
            parts = line.strip().split("\t")
            prot, features = parts[0], [float(f) for f in parts[1:]]
            x[node_dict[prot]] = torch.tensor(features)"""

    # Chargement des attributs d'arêtes
    edge_attr_dict = {}
    with open(VIEW_DIR + edge_attr_file) as f:
        for line in f:
            u, v, val = line.strip().split()
            uid, vid = node_dict[u], node_dict[v]
            edge_attr_dict[(uid, vid)] = float(val)
            edge_attr_dict[(vid, uid)] = float(val)  # facultatif si graphe non orienté

    edge_attr = []
    for i in range(edge_index.size(1)):
        src, tgt = edge_index[0, i].item(), edge_index[1, i].item()
        edge_attr.append(edge_attr_dict.get((src, tgt), 0.0))  # valeur par défaut = 0.0

    edge_attr = torch.tensor(edge_attr, dtype=torch.float).unsqueeze(1)  # shape [E, 1]

    # Normalisation des features de nœuds
    scaler = StandardScaler()
    x = torch.tensor(scaler.fit_transform(x.numpy()), dtype=torch.float)

    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)

    return data, node_dict


class GraphSAGE(torch.nn.Module):
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


# Fonction de perte contrastive InfoNCE
def contrastive_loss(embeddings, edge_index):
    pos_i, pos_j = edge_index  # Arêtes positives (nœuds connectés)
    neg_j = negative_sampling(edge_index, num_nodes=embeddings.size(0))[1]  # Nœuds négatifs

    pos_sim = (embeddings[pos_i] * embeddings[pos_j]).sum(dim=1)
    neg_sim = (embeddings[pos_i] * embeddings[neg_j]).sum(dim=1)
    loss = -torch.log(torch.sigmoid(pos_sim - neg_sim) + 1e-10).mean()  # InfoNCE
    return loss


def train(model, optimizer, x, edge_index, epochs=50):
    for epoch in range(epochs):
        optimizer.zero_grad()
        embeddings = model(x, edge_index)
        loss = contrastive_loss(embeddings, edge_index)
        loss.backward()
        optimizer.step()
        if epoch % 10 == 0:
            print(f"Epoch {epoch}, Loss: {loss.item():.4f}")
            # print(optimizer.state_dict)


def save_embs(model, data, out_file):
    embeddings = model(data.x, data.edge_index)
    node_dict = {}  # Ce dictionnaire doit mapper l'index à l'ID du nœud
    with open(RAWDATA_DIR + "PPInetwork_proteins.txt") as f:
        node_list = f.readlines()
        node_list = [node.split("\n")[0] for node in node_list]
        for index, n in enumerate(node_list):
            node_dict[index] = n
    # print("Emb dim: ", embeddings.shape)
    with open(VIEW_DIR + out_file, "w") as f:
        for i in range(embeddings.size(0)):
            node_name = node_dict[i]  # Récupère le nom du nœud
            embedding_values = " ".join(map(str, embeddings[i].tolist()))  # Convertit les embeddings en string
            f.write(f"{node_name} {embedding_values}\n")

    # print("Embeddings sauvegardés dans 'embeddings.txt'.")


def run_graphsage_pipeline(filename, embed_output_name, hidden_dim=64, out_dim=128, epochs=300, lr=0.01):
    data, _ = load_data(filename)
    model = GraphSAGE(in_channels=data.num_node_features,
                      hidden_channels=hidden_dim,
                      out_channels=out_dim)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    train(model, optimizer, data.x, data.edge_index, epochs=epochs)
    save_embs(model, data, embed_output_name)


if __name__ == '__main__':
    run_graphsage_pipeline("GO-BP_PPInetwork.txt", "GO-BP_PPInetwork_embed.txt")
    run_graphsage_pipeline("GO-CC_PPInetwork.txt", "GO-CC_PPInetwork_embed.txt")
