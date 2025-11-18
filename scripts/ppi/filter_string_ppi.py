import pandas as pd
import networkx as nx

# Load your data
genes_df = pd.read_csv('datasets\\filtered\\consensus_run2\\consensus_genes.csv')
string_df = pd.read_csv('datasets\\ppi\\string_interactions.tsv', sep='\t')

# Your significant genes
significant_genes = set(genes_df['GeneSymbol'])
print(f"Starting with {len(significant_genes)} significant genes")

# STRATEGY 1: Find interactions WHERE AT LEAST ONE gene is in your list
als_network_expanded = string_df[
    string_df['#node1'].isin(significant_genes) | 
    string_df['node2'].isin(significant_genes)
]

print(f"Found {len(als_network_expanded)} interactions involving your genes (at least one side)")

# STRATEGY 2: Use higher confidence threshold to get more reliable connections
als_network_high_conf = als_network_expanded[als_network_expanded['combined_score'] >= 0.7]
print(f"High-confidence interactions (score >= 0.7): {len(als_network_high_conf)}")

def create_expanded_network(interactions_df, seed_genes):
    """Create network including first neighbors of seed genes"""
    G = nx.Graph()
    
    # Add all seed genes as nodes
    for gene in seed_genes:
        G.add_node(gene, type='seed')
    
    # Add interactions and their partners
    for _, row in interactions_df.iterrows():
        G.add_node(row['#node1'], type='seed' if row['#node1'] in seed_genes else 'neighbor')
        G.add_node(row['node2'], type='seed' if row['node2'] in seed_genes else 'neighbor')
        G.add_edge(
            row['#node1'], 
            row['node2'],
            combined_score=row['combined_score'],
            coexpression=row['coexpression']
        )
    
    return G

# Create expanded network
G_expanded = create_expanded_network(als_network_high_conf, significant_genes)

# Network statistics
print(f"\n=== EXPANDED NETWORK ===")
print(f"Total nodes: {G_expanded.number_of_nodes()}")
print(f"Total edges: {G_expanded.number_of_edges()}")
print(f"Seed genes in network: {len([n for n in G_expanded.nodes() if G_expanded.nodes[n]['type'] == 'seed'])}")
print(f"Neighbor genes added: {len([n for n in G_expanded.nodes() if G_expanded.nodes[n]['type'] == 'neighbor'])}")
print(f"Network density: {nx.density(G_expanded):.4f}")

# Analyze the seed genes in the expanded context
def analyze_expanded_network(G, seed_genes, gene_counts):
    """Analyze network focusing on seed genes"""
    degree_centrality = nx.degree_centrality(G)
    betweenness_centrality = nx.betweenness_centrality(G)
    
    results = []
    for gene in seed_genes:
        if gene in G.nodes():
            results.append({
                'GeneSymbol': gene,
                'Type': 'seed',
                'Degree': G.degree(gene),
                'DegreeCentrality': degree_centrality.get(gene, 0),
                'BetweennessCentrality': betweenness_centrality.get(gene, 0),
                'DatasetCount': gene_counts.get(gene, 0),
                'Neighbors': [n for n in G.neighbors(gene) if n not in seed_genes]  # Non-seed neighbors
            })
    
    return pd.DataFrame(results)

# Analyze
gene_count_map = dict(zip(genes_df['GeneSymbol'], genes_df['DatasetCount']))
network_analysis = analyze_expanded_network(G_expanded, significant_genes, gene_count_map)

# Sort by importance
network_analysis['ImportanceScore'] = (
    network_analysis['DegreeCentrality'] + 
    network_analysis['BetweennessCentrality']
)

print("\nTop 10 most connected seed genes:")
top_genes = network_analysis.sort_values('ImportanceScore', ascending=False)
print(top_genes.head(10)[['GeneSymbol', 'Degree', 'ImportanceScore']])

# Check which seed genes are still isolated
connected_seeds = network_analysis[network_analysis['Degree'] > 0]
isolated_seeds = network_analysis[network_analysis['Degree'] == 0]

print(f"\nConnected seed genes: {len(connected_seeds)}")
print(f"Isolated seed genes: {len(isolated_seeds)}")

# Check known ALS genes again
def check_known_als_genes_expanded(G, seed_genes):
    known_als_genes = {
        'SOD1', 'TARDBP', 'FUS', 'C9orf72', 'TBK1', 'OPTN', 'VCP', 
        'UBQLN2', 'SQSTM1', 'PFN1', 'HNRNPA1', 'MATR3'
    }
    
    all_network_genes = set(G.nodes())
    found_als_genes = all_network_genes.intersection(known_als_genes)
    seed_als_genes = seed_genes.intersection(known_als_genes)
    
    print(f"\nKnown ALS genes in entire network: {len(found_als_genes)}/{len(known_als_genes)}")
    print(f"Known ALS genes in your seed list: {len(seed_als_genes)}/{len(known_als_genes)}")
    if found_als_genes:
        print(f"Found in network: {found_als_genes}")
    
    return found_als_genes

known_in_network = check_known_als_genes_expanded(G_expanded, significant_genes)

# Export the expanded network data
# Save all interactions
expanded_interactions = []
for edge in G_expanded.edges(data=True):
    expanded_interactions.append({
        'node1': edge[0],
        'node2': edge[1],
        'node1_type': G_expanded.nodes[edge[0]]['type'],
        'node2_type': G_expanded.nodes[edge[1]]['type'],
        'combined_score': edge[2]['combined_score']
    })

expanded_df = pd.DataFrame(expanded_interactions)
expanded_df.to_csv('als_expanded_network.csv', index=False)
network_analysis.to_csv('als_expanded_analysis.csv', index=False)

print(f"\n=== EXPANDED ANALYSIS COMPLETE ===")
print(f"Now you have a network with {G_expanded.number_of_nodes()} genes and {G_expanded.number_of_edges()} interactions")
print("This gives you a biological context to analyze your significant genes!")