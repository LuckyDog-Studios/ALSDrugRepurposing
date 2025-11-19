import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

def load_existing_network():
    """Load your existing network data"""
    network_df = pd.read_csv('datasets\ppi\\als_expanded_network.csv')
    analysis_df = pd.read_csv('datasets\\ppi\\als_expanded_analysis.csv')
    
    print("=== LOADING EXISTING NETWORK ===")
    print(f"Network edges: {len(network_df)}")
    print(f"Seed genes analyzed: {len(analysis_df)}")
    
    return network_df, analysis_df

def create_enhanced_network(network_df, analysis_df, confidence_threshold=0.4):
    """Create a more complete network using multiple strategies"""
    
    # Get seed genes from analysis
    seed_genes = set(analysis_df['GeneSymbol'])
    print(f"\nSeed genes: {len(seed_genes)}")
    
    # Strategy 1: Use lower confidence threshold to include more interactions
    filtered_network = network_df[network_df['combined_score'] >= confidence_threshold]
    print(f"Edges with confidence >= {confidence_threshold}: {len(filtered_network)}")
    
    # Create initial graph
    G = nx.Graph()
    
    # Add all nodes with their types
    all_nodes = set(filtered_network['node1']).union(set(filtered_network['node2']))
    
    for node in all_nodes:
        node_type = 'seed' if node in seed_genes else 'neighbor'
        G.add_node(node, type=node_type)
    
    # Add edges
    for _, row in filtered_network.iterrows():
        G.add_edge(
            row['node1'], 
            row['node2'],
            combined_score=row['combined_score']
        )
    
    return G, seed_genes

def expand_network_with_ppi(G, seed_genes, ppi_file='datasets\\ppi\\string_interactions.tsv', 
                           confidence_threshold=0.4, max_additions=1000):
    """Expand network by adding relevant PPI interactions"""
    
    try:
        # Load full STRING PPI database
        string_df = pd.read_csv(ppi_file, sep='\t')
        print(f"\nLoaded STRING database: {len(string_df)} interactions")
        
        # Get current network genes
        network_genes = set(G.nodes())
        
        # Find interactions between current network genes that are missing
        missing_interactions = string_df[
            (string_df['#node1'].isin(network_genes)) & 
            (string_df['node2'].isin(network_genes)) &
            (string_df['combined_score'] >= confidence_threshold) &
            ~string_df.apply(lambda x: G.has_edge(x['#node1'], x['node2']), axis=1)
        ]
        
        print(f"Found {len(missing_interactions)} missing interactions between network genes")
        
        # Add missing interactions (limit to prevent explosion)
        added = 0
        for _, row in missing_interactions.head(max_additions).iterrows():
            if not G.has_edge(row['#node1'], row['node2']):
                G.add_edge(
                    row['#node1'], 
                    row['node2'],
                    combined_score=row['combined_score'],
                    source='added_ppi'
                )
                added += 1
        
        print(f"Added {added} missing interactions")
        
        # Add first neighbors of isolated seed genes
        isolated_seeds = [node for node in seed_genes if node in G and G.degree(node) == 0]
        print(f"Isolated seeds to connect: {len(isolated_seeds)}")
        
        for seed in isolated_seeds:
            # Find interactions for isolated seeds
            seed_interactions = string_df[
                ((string_df['#node1'] == seed) | (string_df['node2'] == seed)) &
                (string_df['combined_score'] >= confidence_threshold)
            ]
            
            for _, row in seed_interactions.iterrows():
                partner = row['node2'] if row['#node1'] == seed else row['#node1']
                
                # Add partner node and connection
                if partner not in G:
                    G.add_node(partner, type='neighbor', source='added_for_isolated')
                
                if not G.has_edge(seed, partner):
                    G.add_edge(seed, partner, combined_score=row['combined_score'], source='rescue_edge')
        
    except FileNotFoundError:
        print("STRING PPI file not found, skipping expansion")
    
    return G

def get_largest_connected_component(G, seed_genes):
    """Extract and analyze the largest connected component"""
    
    components = list(nx.connected_components(G))
    component_sizes = [len(c) for c in components]
    
    print(f"\n=== CONNECTED COMPONENTS ===")
    print(f"Number of components: {len(components)}")
    print(f"Component sizes: {component_sizes}")
    
    if not components:
        return G
    
    # Get largest component
    largest_component = max(components, key=len)
    G_lcc = G.subgraph(largest_component).copy()
    
    # Count seed genes in largest component
    seeds_in_lcc = [node for node in G_lcc.nodes() if G_lcc.nodes[node]['type'] == 'seed']
    
    print(f"Largest component: {G_lcc.number_of_nodes()} nodes, {G_lcc.number_of_edges()} edges")
    print(f"Seed genes in largest component: {len(seeds_in_lcc)}/{len(seed_genes)} ({len(seeds_in_lcc)/len(seed_genes)*100:.1f}%)")
    
    return G_lcc

def analyze_enhanced_network(G, seed_genes, original_analysis):
    """Comprehensive analysis of the enhanced network"""
    
    print(f"\n=== ENHANCED NETWORK ANALYSIS ===")
    print(f"Total nodes: {G.number_of_nodes()}")
    print(f"Total edges: {G.number_of_edges()}")
    print(f"Network density: {nx.density(G):.6f}")
    
    # Node type counts
    seed_nodes = [n for n in G.nodes() if G.nodes[n]['type'] == 'seed']
    neighbor_nodes = [n for n in G.nodes() if G.nodes[n]['type'] == 'neighbor']
    
    print(f"Seed nodes: {len(seed_nodes)}")
    print(f"Neighbor nodes: {len(neighbor_nodes)}")
    
    # Connectivity analysis
    connected_seeds = [n for n in seed_nodes if G.degree(n) > 0]
    isolated_seeds = [n for n in seed_nodes if G.degree(n) == 0]
    
    print(f"Connected seed genes: {len(connected_seeds)}")
    print(f"Isolated seed genes: {len(isolated_seeds)}")
    
    # Degree distribution
    degrees = [G.degree(n) for n in G.nodes()]
    print(f"Average degree: {np.mean(degrees):.2f} ± {np.std(degrees):.2f}")
    print(f"Max degree: {max(degrees)}")
    
    # Centrality measures
    degree_centrality = nx.degree_centrality(G)
    betweenness_centrality = nx.betweenness_centrality(G)
    
    # Create enhanced analysis
    enhanced_analysis = []
    for gene in seed_genes:
        if gene in G:
            enhanced_analysis.append({
                'GeneSymbol': gene,
                'Type': 'seed',
                'Degree': G.degree(gene),
                'DegreeCentrality': degree_centrality.get(gene, 0),
                'BetweennessCentrality': betweenness_centrality.get(gene, 0),
                'DatasetCount': original_analysis[original_analysis['GeneSymbol'] == gene]['DatasetCount'].values[0] if gene in original_analysis['GeneSymbol'].values else 0,
                'InLargestComponent': gene in G.nodes(),
                'Neighbors': list(G.neighbors(gene))
            })
        else:
            enhanced_analysis.append({
                'GeneSymbol': gene,
                'Type': 'seed',
                'Degree': 0,
                'DegreeCentrality': 0,
                'BetweennessCentrality': 0,
                'DatasetCount': original_analysis[original_analysis['GeneSymbol'] == gene]['DatasetCount'].values[0] if gene in original_analysis['GeneSymbol'].values else 0,
                'InLargestComponent': False,
                'Neighbors': []
            })
    
    enhanced_df = pd.DataFrame(enhanced_analysis)
    enhanced_df['ImportanceScore'] = (
        enhanced_df['DegreeCentrality'] + 
        enhanced_df['BetweennessCentrality']
    )
    
    return enhanced_df, connected_seeds, isolated_seeds

def identify_key_connectors(G, enhanced_analysis, top_n=20):
    """Identify key connector genes in the network"""
    
    print(f"\n=== KEY CONNECTOR GENES ===")
    
    # Sort by importance
    top_genes = enhanced_analysis.sort_values('ImportanceScore', ascending=False).head(top_n)
    
    print("Top connector seed genes:")
    for _, row in top_genes.head(10).iterrows():
        print(f"  {row['GeneSymbol']}: Degree={row['Degree']}, Importance={row['ImportanceScore']:.4f}")
    
    # Also look at high-degree neighbor genes
    neighbor_genes = [n for n in G.nodes() if G.nodes[n]['type'] == 'neighbor']
    neighbor_degrees = [(n, G.degree(n)) for n in neighbor_genes]
    neighbor_degrees.sort(key=lambda x: x[1], reverse=True)
    
    print(f"\nTop connector neighbor genes:")
    for gene, degree in neighbor_degrees[:10]:
        print(f"  {gene}: Degree={degree}")
    
    return top_genes

def export_enhanced_network(G, enhanced_analysis, filename_prefix='als_enhanced'):
    """Export the enhanced network and analysis"""
    
    # Export network edges
    edges_data = []
    for edge in G.edges(data=True):
        edges_data.append({
            'node1': edge[0],
            'node2': edge[1],
            'node1_type': G.nodes[edge[0]]['type'],
            'node2_type': G.nodes[edge[1]]['type'],
            'combined_score': edge[2].get('combined_score', 0.5),
            'source': edge[2].get('source', 'original')
        })
    
    edges_df = pd.DataFrame(edges_data)
    edges_df.to_csv(f'{filename_prefix}_network.csv', index=False)
    
    # Export analysis
    enhanced_analysis.to_csv(f'{filename_prefix}_analysis.csv', index=False)
    
    # Export node list
    nodes_data = []
    for node in G.nodes(data=True):
        nodes_data.append({
            'GeneSymbol': node[0],
            'Type': node[1]['type'],
            'Degree': G.degree(node[0])
        })
    
    nodes_df = pd.DataFrame(nodes_data)
    nodes_df.to_csv(f'{filename_prefix}_nodes.csv', index=False)
    
    print(f"\n=== EXPORT COMPLETE ===")
    print(f"Network edges: {filename_prefix}_network.csv")
    print(f"Gene analysis: {filename_prefix}_analysis.csv")
    print(f"Node list: {filename_prefix}_nodes.csv")

def main():
    """Main function to create enhanced network"""
    
    # Load existing data
    network_df, analysis_df = load_existing_network()
    
    # Create enhanced network with lower confidence threshold
    print("\n=== CREATING ENHANCED NETWORK ===")
    G, seed_genes = create_enhanced_network(network_df, analysis_df, confidence_threshold=0.4)
    
    print(f"Initial enhanced network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Expand with additional PPI data
    G = expand_network_with_ppi(G, seed_genes, confidence_threshold=0.4)
    
    print(f"After PPI expansion: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Focus on largest connected component
    G_final = get_largest_connected_component(G, seed_genes)
    
    # Analyze the final network
    enhanced_analysis, connected_seeds, isolated_seeds = analyze_enhanced_network(
        G_final, seed_genes, analysis_df
    )
    
    # Identify key connectors
    key_connectors = identify_key_connectors(G_final, enhanced_analysis)
    
    # Compare with original
    original_isolated = len(analysis_df[analysis_df['Degree'] == 0])
    new_isolated = len(isolated_seeds)
    
    print(f"\n=== IMPROVEMENT SUMMARY ===")
    print(f"Original isolated seeds: {original_isolated}/{len(analysis_df)} ({original_isolated/len(analysis_df)*100:.1f}%)")
    print(f"Enhanced isolated seeds: {new_isolated}/{len(seed_genes)} ({new_isolated/len(seed_genes)*100:.1f}%)")
    print(f"Connectivity improvement: {((original_isolated - new_isolated)/original_isolated)*100:.1f}% reduction in isolated seeds")
    
    # Export results
    export_enhanced_network(G_final, enhanced_analysis)
    
    print(f"\n=== ENHANCEMENT COMPLETE ===")
    print(f"Final network: {G_final.number_of_nodes()} genes, {G_final.number_of_edges()} interactions")
    print(f"Connected seed genes: {len(connected_seeds)}")
    print(f"Isolated seed genes: {len(isolated_seeds)}")
    
    if len(isolated_seeds) > 0:
        print(f"Isolated genes: {sorted(isolated_seeds)}")

if __name__ == "__main__":
    main()