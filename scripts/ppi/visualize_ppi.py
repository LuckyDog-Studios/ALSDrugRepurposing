import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

def load_existing_network(analysis_file, network_file):
    """Load existing network analysis and create network graph"""
    print("Loading existing network data...")
    
    # Load your analysis and network files
    analysis_df = pd.read_csv(analysis_file)
    network_df = pd.read_csv(network_file)
    
    # Create network graph from existing data
    G = nx.Graph()
    
    # Add nodes from analysis
    for _, row in analysis_df.iterrows():
        G.add_node(row['GeneSymbol'], type='seed', 
                  degree=row['Degree'],
                  dataset_count=row['DatasetCount'])
    
    # Add edges from network file
    for _, row in network_df.iterrows():
        # Make sure both nodes exist in the graph
        if row['node1'] not in G.nodes():
            G.add_node(row['node1'], type=row['node1_type'])
        if row['node2'] not in G.nodes():
            G.add_node(row['node2'], type=row['node2_type'])
        
        G.add_edge(row['node1'], row['node2'], 
                  combined_score=row['combined_score'])
    
    print(f"Loaded network with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    return G, analysis_df, network_df

def plot_basic_network(G, analysis_df):
    """Create a basic network visualization"""
    print("\nCreating basic network visualization...")
    
    plt.figure(figsize=(14, 12))
    
    # Define node colors by type
    node_colors = []
    node_sizes = []
    for node in G.nodes():
        if G.nodes[node]['type'] == 'seed':
            node_colors.append('red')
            # Larger size for seed genes
            node_sizes.append(G.degree(node) * 100 + 300)
        else:
            node_colors.append('lightblue')
            node_sizes.append(G.degree(node) * 50 + 100)
    
    # Spring layout
    pos = nx.spring_layout(G, k=2, iterations=100, seed=42)
    
    # Draw the network
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, 
                          alpha=0.9, edgecolors='black', linewidths=0.5)
    nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='gray', width=1.5)
    
    # Only label important nodes
    labels = {}
    for node in G.nodes():
        if G.degree(node) > 2 or G.nodes[node]['type'] == 'seed':
            labels[node] = node
    
    nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold')
    
    plt.title("ALS PPI Network\n(Red: Your Significant Genes, Blue: STRING Neighbors)", size=16)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('als_network_basic.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_community_network(G):
    """Visualize with communities colored differently"""
    print("\nCreating community visualization...")
    
    from networkx.algorithms import community
    
    # Detect communities
    communities = community.louvain_communities(G, seed=42)
    
    plt.figure(figsize=(15, 13))
    
    # Create color map for communities
    community_colors = plt.cm.Set3(np.linspace(0, 1, len(communities)))
    node_colors = []
    node_sizes = []
    
    for node in G.nodes():
        for i, comm in enumerate(communities):
            if node in comm:
                node_colors.append(community_colors[i])
                break
        
        # Size by degree and type
        size = G.degree(node) * 80 + 100
        if G.nodes[node]['type'] == 'seed':
            size += 200
        node_sizes.append(size)
    
    pos = nx.spring_layout(G, k=2, iterations=100, seed=42)
    
    # Draw
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, 
                          alpha=0.9, edgecolors='black', linewidths=0.5)
    nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color='gray', width=1.5)
    
    # Only label important nodes
    labels = {}
    for node in G.nodes():
        if G.degree(node) > 3 or G.nodes[node]['type'] == 'seed':
            labels[node] = node
    
    nx.draw_networkx_labels(G, pos, labels, font_size=9, font_weight='bold')
    
    plt.title("ALS Network with Community Detection", size=16)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('als_network_communities.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Found {len(communities)} communities in the network")
    return communities

def create_interactive_network_fixed(G, output_file='als_network_interactive.html'):
    """Create an interactive web-based network visualization (fixed version)"""
    print(f"\nCreating interactive network: {output_file}")
    
    try:
        from pyvis.network import Network
        
        # Create pyvis network
        net = Network(height='800px', width='100%', bgcolor='#ffffff', font_color='black')
        
        # Configure physics for better layout
        net.set_options("""
        var options = {
          "physics": {
            "enabled": true,
            "stabilization": {"iterations": 100},
            "barnesHut": {
              "gravitationalConstant": -80000,
              "springConstant": 0.001,
              "springLength": 200
            }
          }
        }
        """)
        
        # Add nodes
        for node in G.nodes():
            node_type = G.nodes[node]['type']
            degree = G.degree(node)
            
            # Customize node appearance
            if node_type == 'seed':
                color = '#ff4444'  # Red for your genes
                size = degree * 5 + 25
                title = f"Seed Gene: {node}<br>Degree: {degree}<br>Type: Your significant gene"
            else:
                color = '#4444ff'  # Blue for neighbors
                size = degree * 3 + 15
                title = f"Neighbor: {node}<br>Degree: {degree}<br>Type: STRING database neighbor"
            
            net.add_node(node, label=node, color=color, size=size, title=title, font={'size': 20})
        
        # Add edges with weights
        for edge in G.edges(data=True):
            weight = edge[2]['combined_score']
            width = weight * 3
            net.add_edge(edge[0], edge[1], value=weight, width=width, title=f"Score: {weight:.3f}")
        
        # Save as interactive HTML - FIXED: Use write_html instead of show
        net.write_html(output_file)
        print(f"✓ Interactive network saved as {output_file}")
        
    except Exception as e:
        print(f"✗ PyVis interactive network failed: {e}")
        print("Creating alternative visualization instead...")
        create_simple_html_network(G, output_file)

def create_simple_html_network(G, output_file='als_network_simple.html'):
    """Create a simple HTML visualization as fallback"""
    print(f"Creating simple HTML network: {output_file}")
    
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>ALS Network Visualization</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            .node { margin: 10px; padding: 10px; border-radius: 5px; display: inline-block; }
            .seed { background-color: #ff4444; color: white; }
            .neighbor { background-color: #4444ff; color: white; }
            .connection { margin-left: 20px; color: #666; }
        </style>
    </head>
    <body>
        <h1>ALS Protein-Protein Interaction Network</h1>
        <h2>Network Summary</h2>
        <p>Total Nodes: """ + str(G.number_of_nodes()) + """</p>
        <p>Total Edges: """ + str(G.number_of_edges()) + """</p>
        <p>Seed Genes (Red): """ + str(len([n for n in G.nodes() if G.nodes[n]['type'] == 'seed'])) + """</p>
        <p>Neighbor Genes (Blue): """ + str(len([n for n in G.nodes() if G.nodes[n]['type'] == 'neighbor'])) + """</p>
        
        <h2>Connected Components</h2>
    """
    
    # Add connected components
    components = list(nx.connected_components(G))
    for i, component in enumerate(components):
        html_content += f"<h3>Component {i+1} (Size: {len(component)})</h3>"
        for node in component:
            node_type = "seed" if G.nodes[node]['type'] == 'seed' else "neighbor"
            html_content += f'<div class="node {node_type}">{node} (Degree: {G.degree(node)})</div>'
    
    html_content += """
        <h2>Top Connected Genes</h2>
    """
    
    # Add top connected genes
    degrees = [(node, G.degree(node)) for node in G.nodes()]
    top_genes = sorted(degrees, key=lambda x: x[1], reverse=True)[:10]
    
    for node, degree in top_genes:
        node_type = "seed" if G.nodes[node]['type'] == 'seed' else "neighbor"
        html_content += f'<div class="node {node_type}">{node} - {degree} connections</div>'
    
    html_content += """
    </body>
    </html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"✓ Simple HTML network saved as {output_file}")

def plot_pathway_focused(G):
    """Focus on the three main pathways"""
    print("\nCreating pathway-focused visualization...")
    
    # Identify the three main components
    proteasome_nodes = set(nx.node_connected_component(G, 'PSMD7')) if 'PSMD7' in G else set()
    calcium_nodes = set(nx.node_connected_component(G, 'PPP3R1')) if 'PPP3R1' in G else set()
    oxidative_nodes = set(nx.node_connected_component(G, 'CAT')) if 'CAT' in G else set()
    
    # Create subgraph of just these pathways
    pathway_nodes = proteasome_nodes.union(calcium_nodes).union(oxidative_nodes)
    if pathway_nodes:
        pathway_graph = G.subgraph(pathway_nodes)
    else:
        print("No pathway hubs found in network")
        return
    
    plt.figure(figsize=(16, 14))
    
    # Color by pathway
    node_colors = []
    node_sizes = []
    for node in pathway_graph.nodes():
        if node in proteasome_nodes:
            node_colors.append('red')
        elif node in calcium_nodes:
            node_colors.append('blue')
        elif node in oxidative_nodes:
            node_colors.append('green')
        else:
            node_colors.append('gray')
        
        # Size by degree and type
        size = pathway_graph.degree(node) * 120 + 200
        if pathway_graph.nodes[node]['type'] == 'seed':
            size += 300
        node_sizes.append(size)
    
    pos = nx.spring_layout(pathway_graph, k=3, iterations=150, seed=42)
    
    # Draw
    nx.draw_networkx_nodes(pathway_graph, pos, node_color=node_colors, 
                          node_size=node_sizes, alpha=0.9, edgecolors='black', linewidths=1)
    nx.draw_networkx_edges(pathway_graph, pos, alpha=0.5, edge_color='gray', width=2)
    nx.draw_networkx_labels(pathway_graph, pos, font_size=10, font_weight='bold')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', label='Proteasome Pathway (PSMD7)'),
        Patch(facecolor='blue', label='Calcium Signaling (PPP3R1)'),
        Patch(facecolor='green', label='Oxidative Stress (CAT)')
    ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize=12)
    
    plt.title("ALS Core Pathways Network", size=18)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('als_network_pathways.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_degree_distribution(G):
    """Show the connectivity distribution"""
    print("\nCreating degree distribution plot...")
    
    degrees = [G.degree(node) for node in G.nodes()]
    
    plt.figure(figsize=(10, 6))
    plt.hist(degrees, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    plt.xlabel('Number of Connections (Degree)', fontsize=12)
    plt.ylabel('Number of Genes', fontsize=12)
    plt.title('Network Connectivity Distribution', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    # Add statistics
    plt.text(0.65, 0.9, f'Total nodes: {len(degrees)}', transform=plt.gca().transAxes)
    plt.text(0.65, 0.85, f'Mean degree: {np.mean(degrees):.2f}', transform=plt.gca().transAxes)
    plt.text(0.65, 0.8, f'Max degree: {max(degrees)}', transform=plt.gca().transAxes)
    
    plt.tight_layout()
    plt.savefig('als_network_degree_dist.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_summary_report(G, analysis_df):
    """Create a summary report of the network"""
    print("\n=== NETWORK SUMMARY ===")
    
    # Count nodes by type
    seed_nodes = [n for n in G.nodes() if G.nodes[n]['type'] == 'seed']
    neighbor_nodes = [n for n in G.nodes() if G.nodes[n]['type'] == 'neighbor']
    
    print(f"Total nodes: {G.number_of_nodes()}")
    print(f"Seed genes (your data): {len(seed_nodes)}")
    print(f"Neighbor genes (STRING): {len(neighbor_nodes)}")
    print(f"Total edges: {G.number_of_edges()}")
    print(f"Network density: {nx.density(G):.4f}")
    
    # Connected components
    components = list(nx.connected_components(G))
    print(f"Connected components: {len(components)}")
    
    # Top connected genes
    connected_seeds = analysis_df[analysis_df['Degree'] > 0]
    print(f"\nConnected seed genes: {len(connected_seeds)}")
    print("Top connected seed genes:")
    top_genes = connected_seeds.nlargest(5, 'Degree')[['GeneSymbol', 'Degree', 'ImportanceScore']]
    print(top_genes.to_string(index=False))
    
    # Check for known ALS genes
    known_als_genes = {'SOD1', 'TARDBP', 'FUS', 'C9orf72', 'TBK1', 'OPTN', 'VCP'}
    als_in_network = known_als_genes.intersection(set(G.nodes()))
    print(f"\nKnown ALS genes in network: {len(als_in_network)}/{len(known_als_genes)}")
    if als_in_network:
        print(f"Found: {als_in_network}")

def main():
    """Main function to run all visualizations from existing CSV files"""
    # File paths - use your existing files
    analysis_file = 'als_enhanced_analysis.csv'
    network_file = 'als_enhanced_network.csv'
    
    try:
        # Step 1: Load existing network data
        G, analysis_df, network_df = load_existing_network(analysis_file, network_file)
        
        # Step 2: Create summary report
        create_summary_report(G, analysis_df)
        
        # Step 3: Create visualizations
        plot_basic_network(G, analysis_df)
        plot_community_network(G)
        create_interactive_network_fixed(G)  # Use the fixed version
        plot_pathway_focused(G)
        plot_degree_distribution(G)
        
        print("\n=== ALL VISUALIZATIONS COMPLETE ===")
        print("Generated files:")
        print("- als_network_basic.png (basic network)")
        print("- als_network_communities.png (community view)")
        print("- als_network_interactive.html (interactive web version)")
        print("- als_network_pathways.png (pathway focus)")
        print("- als_network_degree_dist.png (connectivity distribution)")
        
    except Exception as e:
        print(f"Error: {e}")
        print("Please make sure your CSV files are in the same directory")

if __name__ == "__main__":
    main()