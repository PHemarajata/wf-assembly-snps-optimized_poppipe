process TREE_GRAFT {
    tag "Tree grafting"
    label "process_medium"

    conda "conda-forge::python=3.9 conda-forge::ete3=3.1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ete3:3.1.3--pyhdfd78af_0' :
        'quay.io/biocontainers/ete3:3.1.3--pyhdfd78af_0' }"

    input:
    path cluster_trees
    path poppipe_output_dir
    path cluster_info

    output:
    path "grafted_tree.nwk", emit: tree
    path "tree_graft_log.txt", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

import sys
import os
import json
import glob
from ete3 import Tree
from pathlib import Path

def midpoint_root(tree):
    \"\"\"Root tree at midpoint\"\"\"
    try:
        R = tree.get_midpoint_outgroup()
        if R:
            tree.set_outgroup(R)
    except:
        pass  # If midpoint rooting fails, continue with original tree
    return tree

def load_cluster_trees():
    \"\"\"Load all cluster trees\"\"\"
    cluster_trees = {}
    tree_files = glob.glob("*.nwk") + glob.glob("*.tree") + glob.glob("*.treefile")
    
    for tree_file in tree_files:
        try:
            tree = Tree(tree_file)
            # Extract cluster ID from filename
            cluster_id = Path(tree_file).stem
            # Remove common prefixes/suffixes
            cluster_id = cluster_id.replace('cluster_', '').replace('_tree', '').replace('_phylogeny', '')
            cluster_trees[cluster_id] = tree
            print(f"Loaded tree for cluster {cluster_id} from {tree_file}")
        except Exception as e:
            print(f"Error loading tree from {tree_file}: {e}")
    
    return cluster_trees

def load_poppipe_nj_trees():
    \"\"\"Load PopPIPE NJ trees for scaling\"\"\"
    nj_trees = {}
    poppipe_dir = "${poppipe_output_dir}"
    
    # Look for NJ trees in strain directories
    strain_dirs = []
    possible_strain_paths = [
        os.path.join(poppipe_dir, 'strain'),
        os.path.join(poppipe_dir, 'strains'),
        os.path.join(poppipe_dir, 'output', 'strain'),
        os.path.join(poppipe_dir, 'output', 'strains')
    ]
    
    for strain_path in possible_strain_paths:
        if os.path.exists(strain_path):
            for cluster_dir in os.listdir(strain_path):
                cluster_path = os.path.join(strain_path, cluster_dir)
                if os.path.isdir(cluster_path):
                    # Look for NJ tree files
                    nj_files = [
                        os.path.join(cluster_path, 'njtree.nwk'),
                        os.path.join(cluster_path, 'nj_tree.nwk'),
                        os.path.join(cluster_path, 'rapidnj.nwk')
                    ]
                    
                    for nj_file in nj_files:
                        if os.path.exists(nj_file):
                            try:
                                nj_tree = Tree(nj_file)
                                nj_trees[cluster_dir] = nj_tree
                                print(f"Loaded NJ tree for cluster {cluster_dir}")
                                break
                            except Exception as e:
                                print(f"Error loading NJ tree from {nj_file}: {e}")
            break
    
    return nj_trees

def load_poppipe_backbone():
    \"\"\"Load PopPIPE backbone tree if available\"\"\"
    backbone_tree = None
    poppipe_dir = "${poppipe_output_dir}"
    
    # Look for backbone tree in PopPIPE output
    possible_backbone_files = [
        os.path.join(poppipe_dir, "full_tree.nwk"),
        os.path.join(poppipe_dir, "backbone.nwk"),
        os.path.join(poppipe_dir, "backbone.tree"),
        os.path.join(poppipe_dir, "output", "full_tree.nwk"),
        os.path.join(poppipe_dir, "strain", "other", "njtree.nwk")
    ]
    
    for backbone_file in possible_backbone_files:
        if os.path.exists(backbone_file):
            try:
                backbone_tree = Tree(backbone_file)
                print(f"Loaded backbone tree from {backbone_file}")
                break
            except Exception as e:
                print(f"Error loading backbone from {backbone_file}: {e}")
    
    return backbone_tree

def scale_tree_branches(ml_tree, nj_tree):
    \"\"\"Scale ML tree branches based on NJ tree total length\"\"\"
    try:
        ml_length = sum(node.dist for node in ml_tree.traverse())
        nj_length = sum(node.dist for node in nj_tree.traverse())
        
        if ml_length > 0:
            scale = nj_length / ml_length
            for node in ml_tree.traverse():
                node.dist *= scale
            print(f"Scaled tree branches by factor {scale:.4f}")
    except Exception as e:
        print(f"Error scaling tree branches: {e}")
    
    return ml_tree

def graft_trees_poppipe_style(backbone_tree, cluster_trees, nj_trees):
    \"\"\"Graft cluster trees using PopPIPE-style algorithm\"\"\"
    
    if backbone_tree is None:
        print("No backbone tree found, creating combined tree from clusters")
        return create_combined_tree(cluster_trees)
    
    print("Using PopPIPE-style tree grafting")
    grafted_tree = backbone_tree.copy()
    
    for cluster_id, ml_tree in cluster_trees.items():
        try:
            # Scale ML tree if corresponding NJ tree exists
            if cluster_id in nj_trees:
                ml_tree = scale_tree_branches(ml_tree.copy(), nj_trees[cluster_id])
            else:
                ml_tree = ml_tree.copy()
            
            # Find a sample from this cluster in the backbone tree
            graft_point = None
            for leaf in ml_tree.get_leaves():
                sample_name = leaf.name
                # Look for this sample in the backbone tree
                backbone_nodes = grafted_tree.search_nodes(name=sample_name)
                if backbone_nodes:
                    graft_point = backbone_nodes[0]
                    break
            
            if graft_point:
                # Set the ML tree root to the sample we found
                ml_tree.set_outgroup(sample_name)
                
                # Replace the graft point with the ML tree
                parent = graft_point.up
                if parent:
                    parent.remove_child(graft_point)
                    parent.add_child(ml_tree, dist=0)
                    print(f"Grafted cluster {cluster_id} tree at sample {sample_name}")
                else:
                    # This is the root
                    grafted_tree = ml_tree
                    print(f"Replaced root with cluster {cluster_id} tree")
            else:
                print(f"Could not find graft point for cluster {cluster_id}")
                
        except Exception as e:
            print(f"Error grafting cluster {cluster_id}: {e}")
            continue
    
    return grafted_tree

def create_combined_tree(cluster_trees):
    \"\"\"Create a combined tree when no backbone is available\"\"\"
    if not cluster_trees:
        raise ValueError("No cluster trees available")
    
    cluster_ids = list(cluster_trees.keys())
    combined_tree = cluster_trees[cluster_ids[0]].copy()
    
    # Add other clusters as sister clades
    for cluster_id in cluster_ids[1:]:
        cluster_tree = cluster_trees[cluster_id].copy()
        
        # Create a new root that includes both trees
        new_root = Tree()
        new_root.add_child(combined_tree)
        new_root.add_child(cluster_tree)
        combined_tree = new_root
    
    print(f"Created combined tree from {len(cluster_trees)} clusters")
    return combined_tree

# Main execution
try:
    print("Starting PopPIPE-style tree grafting process...")
    
    # Load cluster trees (ML trees from SNP analysis)
    cluster_trees = load_cluster_trees()
    print(f"Loaded {len(cluster_trees)} cluster ML trees")
    
    # Load PopPIPE NJ trees for scaling
    nj_trees = load_poppipe_nj_trees()
    print(f"Loaded {len(nj_trees)} PopPIPE NJ trees")
    
    # Load backbone tree
    backbone_tree = load_poppipe_backbone()
    
    # Perform grafting
    if cluster_trees:
        grafted_tree = graft_trees_poppipe_style(backbone_tree, cluster_trees, nj_trees)
        
        # Root the tree at midpoint
        grafted_tree = midpoint_root(grafted_tree)
        
        # Write the grafted tree
        grafted_tree.write(format=5, outfile="grafted_tree.nwk")
        print("Successfully created grafted tree")
        
        # Write log
        with open("tree_graft_log.txt", "w") as log_file:
            log_file.write(f"PopPIPE-style tree grafting completed successfully\\n")
            log_file.write(f"Number of cluster ML trees: {len(cluster_trees)}\\n")
            log_file.write(f"Number of PopPIPE NJ trees: {len(nj_trees)}\\n")
            log_file.write(f"Backbone tree used: {'Yes' if backbone_tree else 'No'}\\n")
            log_file.write(f"Final tree leaves: {len(grafted_tree.get_leaves())}\\n")
            log_file.write(f"Cluster IDs processed: {', '.join(cluster_trees.keys())}\\n")
    
    else:
        raise ValueError("No cluster trees found for grafting")

except Exception as e:
    print(f"Error in tree grafting: {e}")
    
    # Write error log
    with open("tree_graft_log.txt", "w") as log_file:
        log_file.write(f"Tree grafting failed: {e}\\n")
    
    # Create empty tree file to prevent pipeline failure
    with open("grafted_tree.nwk", "w") as f:
        f.write("();\\n")
    
    sys.exit(1)

# Write versions
with open('versions.yml', 'w') as f:
    f.write('"${task.process}":\\n')
    f.write('    python: "3.9"\\n')
    f.write('    ete3: "3.1.3"\\n')
    """
}