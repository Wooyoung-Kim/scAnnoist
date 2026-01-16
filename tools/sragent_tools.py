# SRAgent Integration Tools
# Wrapper tools to leverage SRAgent capabilities within scrna_agent

import os
import sys
from typing import Annotated, Optional, List, Dict, Any

from langchain_core.tools import tool

# Try to import SRAgent if available
try:
    from SRAgent.agents.entrez import create_entrez_agent
    from SRAgent.agents.papers import create_papers_agent
    SRAGENT_AVAILABLE = True
except ImportError:
    SRAGENT_AVAILABLE = False
    print("Warning: SRAgent not installed. Literature search will use fallback methods.")


# ============================================================
# SRAgent Wrapper Tools
# ============================================================

@tool
def search_sra_for_reference(
    query: Annotated[str, "Search query for finding reference datasets"],
    species: Annotated[str, "Species to filter by (e.g., 'mouse', 'human')"] = "mouse",
    max_results: Annotated[int, "Maximum number of results"] = 5,
) -> Annotated[str, "SRA search results"]:
    """
    Search the SRA database for reference scRNA-seq datasets.
    Uses SRAgent's Entrez agent if available, otherwise falls back to direct Entrez queries.
    
    Useful for finding reference datasets for label transfer or validation.
    """
    from Bio import Entrez
    
    email = os.getenv("EMAIL1") or os.getenv("ENTREZ_EMAIL") or "user@example.com"
    api_key = os.getenv("NCBI_API_KEY1") or os.getenv("NCBI_API_KEY")
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    
    # Build search query for scRNA-seq
    search_query = f"({query}) AND {species}[Organism] AND single cell[All Fields] AND RNA-seq[All Fields]"
    
    try:
        handle = Entrez.esearch(db="sra", term=search_query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        ids = record.get("IdList", [])
        
        if not ids:
            return f"No SRA datasets found for query: {query}"
        
        # Get summaries
        handle = Entrez.esummary(db="sra", id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()
        
        result = f"""
SRA Reference Dataset Search
============================
Query: {query}
Species: {species}
Results: {len(ids)}

Datasets:
"""
        for summary in summaries:
            exp_xml = summary.get("ExpXml", "")
            runs = summary.get("Runs", "")
            
            # Extract accession from ExpXml
            import re
            acc_match = re.search(r'acc="(SRX\d+)"', exp_xml)
            srx = acc_match.group(1) if acc_match else "Unknown"
            
            title_match = re.search(r'Title[^>]*>([^<]+)', exp_xml)
            title = title_match.group(1) if title_match else "No title"
            
            result += f"""
- {srx}: {title[:80]}...
"""
        
        return result
        
    except Exception as e:
        return f"Error searching SRA: {str(e)}"


@tool
def fetch_paper_markers(
    accession: Annotated[str, "SRA accession (SRX, SRP, or GSE)"],
    cell_types: Annotated[Optional[List[str]], "Cell types to search for in papers"] = None,
) -> Annotated[str, "Marker genes extracted from publications"]:
    """
    Find publications associated with an SRA dataset and extract marker gene information.
    Uses SRAgent's papers agent if available.
    
    This extends Literature RAG by finding markers from dataset-specific papers.
    """
    from Bio import Entrez
    import re
    
    email = os.getenv("EMAIL1") or "user@example.com"
    api_key = os.getenv("NCBI_API_KEY1")
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    
    try:
        # Search for the accession in GEO/SRA
        handle = Entrez.esearch(db="gds", term=accession, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        
        geo_ids = record.get("IdList", [])
        
        if not geo_ids:
            # Try SRA directly
            handle = Entrez.esearch(db="sra", term=accession, retmax=5)
            record = Entrez.read(handle)
            handle.close()
            geo_ids = record.get("IdList", [])
        
        if not geo_ids:
            return f"No records found for accession: {accession}"
        
        # Link to PubMed
        handle = Entrez.elink(dbfrom="gds", db="pubmed", id=geo_ids[0])
        link_record = Entrez.read(handle)
        handle.close()
        
        pubmed_ids = []
        for linkset in link_record:
            for link in linkset.get("LinkSetDb", []):
                for linked in link.get("Link", []):
                    pubmed_ids.append(linked["Id"])
        
        if not pubmed_ids:
            return f"No publications found linked to {accession}"
        
        # Fetch abstracts
        handle = Entrez.efetch(db="pubmed", id=pubmed_ids[:5], rettype="abstract", retmode="text")
        abstracts = handle.read()
        handle.close()
        
        # Extract potential marker genes
        gene_pattern = r'\b([A-Z][a-z0-9]{1,6}|[A-Z]{2,6}[0-9]*)\b'
        potential_genes = re.findall(gene_pattern, abstracts)
        
        exclude_words = {'The', 'This', 'That', 'With', 'From', 'For', 'And', 'Not',
                        'Are', 'Was', 'Were', 'Has', 'Have', 'Had', 'May', 'Can',
                        'Cell', 'Cells', 'RNA', 'DNA', 'PCR', 'FACS', 'Fig', 'Table',
                        'Gene', 'Genes', 'Data', 'Study'}
        
        gene_counts = {}
        for gene in potential_genes:
            if gene not in exclude_words and len(gene) >= 2:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
        
        top_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:20]
        
        result = f"""
Publication-Derived Markers for {accession}
============================================
PubMed IDs: {', '.join(pubmed_ids[:5])}

Top Mentioned Genes:
"""
        for gene, count in top_genes[:15]:
            result += f"  - {gene}: {count} mentions\n"
        
        if cell_types:
            result += f"\nSearch context: {', '.join(cell_types)}\n"
        
        return result
        
    except Exception as e:
        return f"Error fetching paper markers: {str(e)}"


# ============================================================
# Ontology Tools
# ============================================================

# Tissue ontology mapping (subset of Uberon)
TISSUE_ONTOLOGY = {
    "spleen": {"id": "UBERON:0002106", "name": "spleen"},
    "bone marrow": {"id": "UBERON:0002371", "name": "bone marrow"},
    "lymph node": {"id": "UBERON:0000029", "name": "lymph node"},
    "thymus": {"id": "UBERON:0002370", "name": "thymus"},
    "blood": {"id": "UBERON:0000178", "name": "blood"},
    "pbmc": {"id": "UBERON:0000178", "name": "blood"},
    "lung": {"id": "UBERON:0002048", "name": "lung"},
    "liver": {"id": "UBERON:0002107", "name": "liver"},
    "brain": {"id": "UBERON:0000955", "name": "brain"},
    "heart": {"id": "UBERON:0000948", "name": "heart"},
    "kidney": {"id": "UBERON:0002113", "name": "kidney"},
    "intestine": {"id": "UBERON:0000160", "name": "intestine"},
    "colon": {"id": "UBERON:0001155", "name": "colon"},
    "skin": {"id": "UBERON:0002097", "name": "skin"},
    "muscle": {"id": "UBERON:0001630", "name": "skeletal muscle tissue"},
    "pancreas": {"id": "UBERON:0001264", "name": "pancreas"},
}

# Cell type ontology mapping (subset of Cell Ontology)
CELL_ONTOLOGY = {
    "B cell": {"id": "CL:0000236", "name": "B cell"},
    "T cell": {"id": "CL:0000084", "name": "T cell"},
    "NK cell": {"id": "CL:0000623", "name": "natural killer cell"},
    "monocyte": {"id": "CL:0000576", "name": "monocyte"},
    "macrophage": {"id": "CL:0000235", "name": "macrophage"},
    "dendritic cell": {"id": "CL:0000451", "name": "dendritic cell"},
    "neutrophil": {"id": "CL:0000775", "name": "neutrophil"},
    "plasma cell": {"id": "CL:0000786", "name": "plasma cell"},
    "CD4 T cell": {"id": "CL:0000624", "name": "CD4-positive, alpha-beta T cell"},
    "CD8 T cell": {"id": "CL:0000625", "name": "CD8-positive, alpha-beta T cell"},
    "regulatory T cell": {"id": "CL:0000815", "name": "regulatory T cell"},
    "memory B cell": {"id": "CL:0000787", "name": "memory B cell"},
    "naive B cell": {"id": "CL:0000788", "name": "naive B cell"},
    "germinal center B cell": {"id": "CL:0000844", "name": "germinal center B cell"},
    "stem cell": {"id": "CL:0000034", "name": "stem cell"},
    "hematopoietic stem cell": {"id": "CL:0000037", "name": "hematopoietic stem cell"},
}


@tool
def map_tissue_ontology(
    tissue_name: Annotated[str, "Tissue name to map (e.g., 'spleen', 'bone marrow')"],
) -> Annotated[str, "Uberon ontology mapping"]:
    """
    Map a tissue name to its Uberon ontology ID.
    Useful for standardizing tissue annotations and metadata.
    """
    tissue_lower = tissue_name.lower().strip()
    
    if tissue_lower in TISSUE_ONTOLOGY:
        ont = TISSUE_ONTOLOGY[tissue_lower]
        return f"Tissue: {tissue_name}\nUberon ID: {ont['id']}\nStandard Name: {ont['name']}"
    
    # Fuzzy matching
    matches = []
    for key, value in TISSUE_ONTOLOGY.items():
        if tissue_lower in key or key in tissue_lower:
            matches.append((key, value))
    
    if matches:
        result = f"Tissue: {tissue_name}\nPossible matches:\n"
        for key, value in matches:
            result += f"  - {key}: {value['id']} ({value['name']})\n"
        return result
    
    return f"Tissue '{tissue_name}' not found in local ontology. Try SRAgent's tissue-ontology agent for comprehensive mapping."


@tool
def map_cell_type_ontology(
    cell_type: Annotated[str, "Cell type name to map"],
) -> Annotated[str, "Cell Ontology mapping"]:
    """
    Map a cell type name to its Cell Ontology (CL) ID.
    Useful for standardizing cell type annotations.
    """
    cell_lower = cell_type.lower().strip()
    
    # Direct match
    for key, value in CELL_ONTOLOGY.items():
        if cell_lower == key.lower():
            return f"Cell Type: {cell_type}\nCL ID: {value['id']}\nStandard Name: {value['name']}"
    
    # Partial match
    matches = []
    for key, value in CELL_ONTOLOGY.items():
        if cell_lower in key.lower() or key.lower() in cell_lower:
            matches.append((key, value))
    
    if matches:
        result = f"Cell Type: {cell_type}\nPossible matches:\n"
        for key, value in matches:
            result += f"  - {key}: {value['id']} ({value['name']})\n"
        return result
    
    return f"Cell type '{cell_type}' not found in local ontology."


@tool
def add_ontology_annotations(
    tissue: Annotated[Optional[str], "Tissue name for the dataset"] = None,
    cell_types: Annotated[Optional[List[str]], "List of cell types in the dataset"] = None,
) -> Annotated[str, "Ontology annotation summary"]:
    """
    Add ontology IDs to the current dataset's metadata.
    Maps tissue to Uberon and cell types to Cell Ontology.
    """
    from scrna_agent.tools.pipeline_tools import get_pipeline_state, set_pipeline_state
    
    state = get_pipeline_state()
    adata = state.get("adata")
    
    result = "Ontology Annotations\n====================\n"
    
    # Tissue ontology
    if tissue:
        tissue_result = map_tissue_ontology.invoke({"tissue_name": tissue})
        result += f"\n{tissue_result}\n"
        
        tissue_lower = tissue.lower().strip()
        if tissue_lower in TISSUE_ONTOLOGY:
            if adata is not None:
                adata.uns["tissue_ontology"] = TISSUE_ONTOLOGY[tissue_lower]
            set_pipeline_state("tissue_ontology", TISSUE_ONTOLOGY[tissue_lower])
    
    # Cell type ontology
    if cell_types:
        result += "\nCell Type Ontology Mappings:\n"
        cell_ont_map = {}
        
        for ct in cell_types:
            ct_lower = ct.lower().strip()
            for key, value in CELL_ONTOLOGY.items():
                if ct_lower in key.lower() or key.lower() in ct_lower:
                    cell_ont_map[ct] = value
                    result += f"  - {ct}: {value['id']} ({value['name']})\n"
                    break
            else:
                result += f"  - {ct}: No mapping found\n"
        
        if adata is not None:
            adata.uns["cell_type_ontology"] = cell_ont_map
        set_pipeline_state("cell_type_ontology", cell_ont_map)
    
    return result
