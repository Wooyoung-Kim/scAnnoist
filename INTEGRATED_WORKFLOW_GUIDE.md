# Integrated scRNA-seq Analysis Workflow with AI Skills

## 🎯 Overview

이 워크플로우는 최신 AI skills를 통합하여 scRNA-seq 분석을 완전히 자동화합니다:

- **Multi-Agent Annotation**: 4명의 전문 agents가 토론하여 cell type 결정
- **Scientific Visualization**: `scientific-visualization` skill 사용
- **Professional PowerPoint**: `pptx` + `scientific-slides` skills 사용

---

## 🚀 Quick Start

```bash
# 전체 워크플로우 실행
python master_workflow.py
```

이 스크립트가 단계별로 안내합니다.

---

## 📋 Complete Workflow

### **Phase 1: Basic Analysis** (자동)

```bash
python JMV_ALI_analysis.py
```

**출력:**
- Processed h5ad file
- Basic figures
- Marker genes (CSV)
- Basic PPT

**소요 시간:** ~10분

---

### **Phase 2: Multi-Agent Annotation** (Interactive)

#### Step 2.1: 준비

```bash
python multi_agent_annotation.py
```

**생성:**
- Agent prompts (4개)
- Cluster information (JSON)
- Workflow instructions

#### Step 2.2: Agent 실행 (Claude Code에서)

**Agent 1: Cell Type Expert**
```
I need a general-purpose agent to analyze clusters and propose cell types.

Read: /mnt2/kwy/scrna_agent/JMV_ALI_results/agent_prompts/cell_type_expert_prompt.txt

Provide detailed annotations in JSON format for all clusters.
```

**Agent 2: Computational Biologist**
```
I need a general-purpose agent to evaluate cluster quality.

Read: /mnt2/kwy/scrna_agent/JMV_ALI_results/agent_prompts/computational_biologist_prompt.txt

Provide quantitative assessment in JSON format.
```

**Agent 3: Literature Curator**
```
I need a general-purpose agent to validate with literature.

Read: /mnt2/kwy/scrna_agent/JMV_ALI_results/agent_prompts/literature_curator_prompt.txt

Provide literature support in JSON format.
```

**Agent 4: Critical Evaluator**
```
I need a general-purpose agent to synthesize all opinions.

First read the outputs from agents 1-3, then read:
/mnt2/kwy/scrna_agent/JMV_ALI_results/agent_prompts/critical_evaluator_prompt.txt

Provide final consensus annotations in JSON format.
```

**소요 시간:** ~30-40분 (agents 병렬 실행 가능)

---

### **Phase 3: Scientific Visualizations** (Skill 사용)

#### Step 3.1: 준비

```bash
python generate_scientific_visualizations.py
```

**생성:**
- visualization_specifications.json
- VISUALIZATION_INSTRUCTIONS.md

#### Step 3.2: Skill 실행 (Claude Code에서)

```bash
/scientific-visualization
```

**제공할 정보:**
- Specifications: `/mnt2/kwy/scrna_agent/JMV_ALI_results/visualization_specifications.json`
- Data: `/mnt2/kwy/scrna_agent/JMV_ALI_results/JMV_ALI_harmony_integrated.h5ad`
- Output: `/mnt2/kwy/scrna_agent/JMV_ALI_results/figures/`

**생성되는 8개 Publication-Quality Figures:**
1. QC Metrics Panel (multi-panel)
2. UMAP Comprehensive (4-panel)
3. Marker Gene Heatmap
4. Cell Type Proportions
5. Marker Gene Dotplot
6. Volcano Plot
7. PCA Biplot
8. Gene Expression Violins

**소요 시간:** ~15-20분

---

### **Phase 4: Scientific PowerPoint** (Skills 사용)

#### Step 4.1: 준비

```bash
python generate_scientific_ppt.py
```

**생성:**
- presentation_outline.json
- ppt_generation_instructions.txt

#### Step 4.2: Skill 실행 (Claude Code에서)

**Option A: scientific-slides skill (권장)**
```bash
/scientific-slides
```

**Option B: pptx skill**
```bash
/pptx
```

**제공할 정보:**
- Outline: `/mnt2/kwy/scrna_agent/JMV_ALI_results/presentation_outline.json`
- Figures: `/mnt2/kwy/scrna_agent/JMV_ALI_results/figures/`
- Instructions: `/mnt2/kwy/scrna_agent/JMV_ALI_results/ppt_generation_instructions.txt`

**생성되는 내용:**
- 25+ professional slides
- All figures integrated
- Multi-agent annotation summary
- Scientific styling
- Publication-ready

**소요 시간:** ~10-15분

---

## 📊 Output Structure

```
JMV_ALI_results/
├── 📊 Data Files
│   ├── JMV_ALI_harmony_integrated.h5ad (7.3 GB)
│   ├── JMV_ALI_harmony_integrated_annotated.h5ad (7.5 GB)
│   └── metadata.csv
│
├── 📄 Reports & Documentation
│   ├── ANNOTATION_REPORT.md (multi-agent consensus)
│   ├── MULTI_AGENT_WORKFLOW.md
│   ├── VISUALIZATION_INSTRUCTIONS.md
│   └── JMV_ALI_ANALYSIS_SUMMARY.md
│
├── 🤖 Agent System
│   ├── cluster_information_for_agents.json
│   ├── agent_prompts/
│   │   ├── cell_type_expert_prompt.txt
│   │   ├── computational_biologist_prompt.txt
│   │   ├── literature_curator_prompt.txt
│   │   └── critical_evaluator_prompt.txt
│   └── agent_outputs/
│       ├── cell_type_expert_output.json
│       ├── computational_biologist_output.json
│       ├── literature_curator_output.json
│       └── critical_evaluator_output.json
│
├── 📈 Visualizations
│   ├── visualization_specifications.json
│   └── figures/
│       ├── qc_metrics_publication.png (300 DPI)
│       ├── umap_comprehensive_publication.png
│       ├── marker_genes_heatmap_publication.png
│       ├── cell_type_proportions_publication.png
│       ├── marker_dotplot_publication.png
│       ├── volcano_plot_publication.png
│       ├── pca_biplot_publication.png
│       └── gene_expression_violins_publication.png
│
├── 🎤 Presentation
│   ├── presentation_outline.json
│   ├── ppt_generation_instructions.txt
│   └── JMV_ALI_Scientific_Presentation.pptx (30+ MB)
│
└── 📋 Marker Genes
    ├── markers_cluster_0.csv
    ├── markers_cluster_1.csv
    └── ... (21 files)
```

---

## 🎯 Key Features

### 1. Multi-Agent Cell Type Annotation ⭐
- **4 Specialized Agents:**
  - Cell Type Expert: Biological interpretation
  - Computational Biologist: Statistical validation
  - Literature Curator: Published evidence
  - Critical Evaluator: Consensus synthesis

- **Benefits:**
  - Multiple perspectives reduce bias
  - Cross-validation increases confidence
  - Complete rationale documentation
  - Transparent decision-making

### 2. Scientific Visualizations (Skill-Generated) ⭐
- **Publication Quality:**
  - 300 DPI resolution
  - Colorblind-safe palettes
  - Professional layouts
  - Vector-compatible

- **8 Comprehensive Figures:**
  - Multi-panel QC metrics
  - Comprehensive UMAP views
  - Hierarchical heatmaps
  - Statistical comparisons
  - Gene expression patterns

### 3. Professional PowerPoint (Skill-Generated) ⭐
- **Scientific Presentation:**
  - 25+ slides
  - Consistent styling
  - Figure integration
  - Agent discussion summaries
  - Ready for conferences

---

## 💡 Workflow Comparison

### Traditional Approach
```
1. Manual analysis (days)
2. Individual cell type annotation (hours)
3. Basic plots in Python (hours)
4. Manual PowerPoint creation (hours)
Total: Several days
Quality: Variable
```

### Integrated AI Workflow
```
1. Automated analysis (10 min)
2. Multi-agent annotation (40 min)
3. Skill-generated visualizations (20 min)
4. Skill-generated PowerPoint (15 min)
Total: ~90 minutes
Quality: Publication-ready
```

---

## 🔧 Customization

### Modify Agent Prompts
Edit files in `agent_prompts/` to:
- Change tissue-specific knowledge
- Adjust annotation criteria
- Add new expert perspectives
- Modify decision logic

### Customize Visualizations
Edit `visualization_specifications.json` to:
- Add new figure types
- Change color schemes
- Adjust layouts
- Modify statistical tests

### Tailor PowerPoint
Edit `presentation_outline.json` to:
- Add/remove slides
- Change section structure
- Modify content focus
- Adjust styling

---

## 📖 Usage Examples

### Example 1: Quick Analysis
```bash
# Just the basics
python JMV_ALI_analysis.py
```

### Example 2: With Annotation
```bash
# Analysis + Multi-agent annotation
python JMV_ALI_analysis.py
python multi_agent_annotation.py
# Then run agents in Claude Code
```

### Example 3: Full Pipeline
```bash
# Everything with skills
python master_workflow.py
# Follow interactive prompts
```

### Example 4: Custom Workflow
```python
# In your own script
from multi_agent_annotation import prepare_cluster_information
from generate_scientific_visualizations import create_visualization_specifications
from generate_scientific_ppt import generate_ppt_with_skills

# Your custom workflow
```

---

## 🎓 Best Practices

### 1. Agent Discussion
- Run agents 1-3 in parallel
- Review each output before final evaluation
- Iterate if disagreements are unclear
- Document unusual findings

### 2. Visualization
- Check specifications before skill execution
- Validate color schemes for colorblind accessibility
- Ensure all labels are readable
- Test figure sizes for publication

### 3. PowerPoint
- Review outline before generation
- Check figure quality in slides
- Verify statistical claims
- Practice presentation timing

### 4. Documentation
- Keep all agent outputs
- Document any manual overrides
- Note biological interpretations
- Save methodology details

---

## 🚨 Troubleshooting

### Agent Discussion Issues
**Problem:** Agents disagree strongly
**Solution:**
- Check marker gene expression levels
- Review QC metrics for cluster
- Consult additional literature
- Consider sub-clustering

### Visualization Issues
**Problem:** Figures look wrong
**Solution:**
- Verify data file is correct
- Check specification JSON syntax
- Ensure all required columns exist
- Validate color mappings

### PowerPoint Issues
**Problem:** Skill doesn't integrate figures
**Solution:**
- Check file paths in outline
- Verify figures exist and are readable
- Ensure consistent naming
- Regenerate outline if needed

---

## 📚 Additional Resources

### Documentation
- `README.md` - Quick start
- `ANALYSIS_GUIDE.md` - Detailed methods
- `ANNOTATION_REPORT.md` - Agent consensus (auto-generated)
- `VISUALIZATION_INSTRUCTIONS.md` - Skill usage

### Scripts
- `master_workflow.py` - Complete pipeline
- `JMV_ALI_analysis.py` - Basic analysis
- `multi_agent_annotation.py` - Agent system
- `generate_scientific_visualizations.py` - Viz specs
- `generate_scientific_ppt.py` - PPT specs

### Skills Used
- `/pptx` - PowerPoint generation
- `/scientific-slides` - Scientific presentations
- `/scientific-visualization` - Publication figures

---

## ⏱️ Timeline Summary

| Phase | Component | Time | Automation |
|-------|-----------|------|------------|
| 1 | Basic Analysis | 10 min | 100% |
| 2a | Agent Preparation | 5 min | 100% |
| 2b | Agent Discussion | 30 min | 90% (need to run) |
| 3a | Viz Preparation | 5 min | 100% |
| 3b | Viz Generation | 15 min | 95% (skill) |
| 4a | PPT Preparation | 5 min | 100% |
| 4b | PPT Generation | 10 min | 95% (skill) |
| **Total** | | **80 min** | **~95%** |

---

## 🎉 Success Criteria

Workflow is complete when you have:
- [ ] Processed h5ad file with annotations
- [ ] ANNOTATION_REPORT.md with multi-agent consensus
- [ ] 8 publication-quality figures (300 DPI)
- [ ] Professional PowerPoint presentation (25+ slides)
- [ ] All agent outputs documented
- [ ] Marker gene tables
- [ ] Comprehensive documentation

---

## 🔮 Future Enhancements

Planned additions:
- Automated agent compilation script
- More visualization types
- Interactive dashboard generation
- Manuscript figure panel generation
- Supplementary table auto-formatting

---

**Last Updated:** 2026-01-15
**Version:** 2.0 - Integrated AI Skills
**Maintained by:** scRNA-seq Analysis Pipeline Team
