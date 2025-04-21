import torch
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv
from torch_geometric.utils import negative_sampling
from torch_geometric.data import Data
from sklearn.preprocessing import StandardScaler

RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/PPInetwork_topology_view/"


def load_data():
    node_dict = {}
    with open(RAWDATA_DIR + "PPInetwork_proteins.txt") as f:
        node_list = f.readlines()
        node_list = [node.split("\n")[0] for node in node_list]
        index = 0
        for n in node_list:
            node_dict[n] = index
            index += 1
    # print("num nodes: ", len(node_dict))
    with open(RAWDATA_DIR + "PPInetwork_edgelist.csv") as f:
        edge_list = f.readlines()
        edge_list.pop(0)
        edge_list = [inter.split(",") for inter in edge_list]
        edge_list = [(u, v.split("\n")[0]) for u, v in edge_list]
        edge_list = [[node_dict[u], node_dict[v]] for u, v in edge_list]
        edge_list = torch.tensor(edge_list)
    # print("num edges: ", len(edge_list))

    x = torch.zeros(size=(len(node_dict), 10))
    with open(VIEW_DIR + "PPInetwork_topology.txt") as f:
        lines = f.readlines()
        for tab in lines[1:]:
            tab = tab.split("\t")
            tab[-1] = tab[-1].split("\n")[0]
            prot = tab[0]
            features = tab[1:]
            features = [float(i) for i in features]
            x[node_dict[prot]] = torch.tensor(features)

    scaler = StandardScaler()
    data = Data(x=x, edge_index=edge_list.t().contiguous())
    data.x = torch.tensor(scaler.fit_transform(data.x.numpy()), dtype=torch.float)

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


def save_embs(model, data):
    embeddings = model(data.x, data.edge_index)
    node_dict = {}  # Ce dictionnaire doit mapper l'index à l'ID du nœud
    with open(RAWDATA_DIR + "PPInetwork_proteins.txt") as f:
        node_list = f.readlines()
        node_list = [node.split("\n")[0] for node in node_list]
        for index, n in enumerate(node_list):
            node_dict[index] = n
    # print("Emb dim: ", embeddings.shape)
    with open(VIEW_DIR + "SAGE_PPInetwork.txt", "w") as f:
        for i in range(embeddings.size(0)):
            node_name = node_dict[i]  # Récupère le nom du nœud
            embedding_values = " ".join(map(str, embeddings[i].tolist()))  # Convertit les embeddings en string
            f.write(f"{node_name} {embedding_values}\n")

    # print("Embeddings sauvegardés dans 'embeddings.txt'.")


if __name__ == '__main__':
    data, node_dict = load_data()
    model = GraphSAGE(in_channels=data.num_node_features, hidden_channels=64, out_channels=128)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    train(model, optimizer, data.x, data.edge_index, epochs=300)
    save_embs(model, data)
