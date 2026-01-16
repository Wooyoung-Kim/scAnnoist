# PPT Presentation Generation Guide

## Overview

scRNA-agent can automatically generate professional PowerPoint presentations from your analysis results. The system uses Claude Code's `scientific-slides` skill powered by Nano Banana Pro AI to create visually stunning slides.

## Features

✅ **Automatic slide generation** from analysis results
✅ **Professional design** with scientific aesthetics
✅ **Visual consistency** across all slides
✅ **Key metrics and plots** automatically included
✅ **PDF output** ready for presenting

---

## Quick Start

### Method 1: CLI Command (Easiest)

```bash
# Run annotation first
scrna-agent annotate \
  --file data.h5ad \
  --tissue-map tissue_map.csv \
  --output annotated.h5ad

# Generate presentation from results
scrna-agent report \
  --output-dir ./celltypist_outputs \
  --title "My scRNA-seq Study"
```

### Method 2: Python API

```python
from scrna_agent.tools.report_generator import generate_report

# Generate presentation
pdf_path = generate_report(
    output_dir="./celltypist_outputs",
    title="scRNA-seq Analysis Results"
)

print(f"Presentation: {pdf_path}")
```

---

## Complete Workflow

### Step 1: Run Analysis

```bash
# Create tissue mapping
cat > tissue_map.csv << EOF
sample,tissue
sample1,Spleen
sample2,Blood
sample3,Liver
EOF

# Run annotation
scrna-agent annotate \
  --file data.h5ad \
  --tissue-map tissue_map.csv \
  --sample-col sample \
  --output annotated.h5ad
```

**Output files created:**
- `celltypist_predictions.csv`
- `tissue_model_usage.csv`
- `run_metadata.json`

### Step 2: Generate Presentation

```bash
# Generate PPT from results
scrna-agent report \
  --output-dir ./celltypist_outputs \
  --title "Tissue-Aware Cell Type Annotation"
```

**Output:**
- `presentation.pdf` - Complete slide deck
- `slides/` - Individual slide images

---

## Presentation Contents

The auto-generated presentation includes:

### Slide 1: Title
- Study title
- Subtitle: "Tissue-Aware Cell Type Annotation"
- Author attribution

### Slide 2: Overview
- Total cells analyzed
- Number of tissues
- Annotation method

### Slide 3: Tissue Distribution
- Bar chart of cells per tissue
- Visual representation

### Slide 4: Model Selection Strategy
- Decision tree/flowchart
- Priority hierarchy explained

### Slide 5: Models Used
- Which model for each tissue
- Selection rules applied

### Slide 6: Cell Type Distribution
- Top 10 cell types identified
- Pie/donut chart visualization

### Slide 7: Confidence Scores
- Mean/median confidence
- Low confidence cell counts
- Gauge charts

### Slide 8: Summary
- Key takeaways
- Next steps

---

## Python API Examples

### Example 1: Basic Usage

```python
from scrna_agent.tools.report_generator import generate_report

# Simple one-liner
pdf = generate_report("./celltypist_outputs", "My Study")
```

### Example 2: Custom Control

```python
from scrna_agent.tools.report_generator import PresentationGenerator

# Create generator
gen = PresentationGenerator("./celltypist_outputs")

# Load results
results = gen.load_analysis_results()

# Create custom slide plan
slides = gen.create_slide_plan(results, "Custom Title")

# Generate
pdf = gen.generate_presentation("Custom Title")
```

### Example 3: Complete Pipeline

```python
from scrna_agent.tools.annotation_tools import *
from scrna_agent.tools.report_generator import generate_report

# 1. Annotate
load_adata_for_annotation("data.h5ad")
add_tissue_metadata(tissue_map="map.csv")
celltypist_annotate_tissue_aware(
    output_dir="./outputs"
)

# 2. Generate presentation
pdf = generate_report("./outputs", "Results")
```

---

## Customization

### Custom Slide Plan

You can customize which slides to include:

```python
from scrna_agent.tools.report_generator import PresentationGenerator

gen = PresentationGenerator("./outputs")
results = gen.load_analysis_results()

# Create custom slides
custom_slides = [
    {
        "type": "title",
        "title": "My Custom Title",
        "prompt": "Create a professional title slide..."
    },
    {
        "type": "content",
        "title": "Custom Content",
        "prompt": "Create a slide showing..."
    }
]

# Generate each slide
for i, slide in enumerate(custom_slides, 1):
    gen.generate_slide(slide, i)

# Combine to PDF
gen.combine_slides_to_pdf("custom.pdf")
```

### Design Customization

Slides use consistent design with:
- **Color scheme**: Dark blue background, white text, gold accents
- **Typography**: Bold sans-serif titles, clean body text
- **Visual style**: Minimal, professional, modern
- **Layout**: Generous white space, left-aligned content

To customize, edit prompts in `report_generator.py`.

---

## Requirements

### Software Dependencies

```bash
# Install required packages
pip install pillow reportlab

# Ensure scientific-slides skill is available
ls skills/scientific-slides/
```

### Input Files Required

The presentation generator expects these files in `output_dir`:

1. **celltypist_predictions.csv**
   - Per-cell predictions
   - Must have: cell_id, sample, tissue, label, confidence

2. **tissue_model_usage.csv**
   - Models used per tissue
   - Must have: tissue, model_used, selection_rule, n_cells

3. **run_metadata.json**
   - Overall statistics
   - Must have: total_cells, total_tissues, tissue_model_mapping

### Nano Banana Pro

The slide generation uses Nano Banana Pro AI. Ensure it's accessible:

```bash
# Test Nano Banana Pro
python skills/scientific-slides/scripts/generate_slide_image_ai.py \
  --prompt "Test slide" \
  --output test.png
```

---

## Troubleshooting

### Issue: "No analysis results found"

**Solution**: Run annotation first and specify correct output directory

```bash
# Check output directory contains required files
ls celltypist_outputs/
# Should show: celltypist_predictions.csv, tissue_model_usage.csv, run_metadata.json
```

### Issue: "Slide generation timeout"

**Solution**: Increase timeout or simplify prompts

```python
# In report_generator.py, line ~XXX
result = subprocess.run(cmd, timeout=300)  # Increase from 120 to 300
```

### Issue: "PDF creation failed"

**Solution**: Check all slides were generated

```bash
ls celltypist_outputs/slides/
# Should show: slide_01.png, slide_02.png, etc.
```

### Issue: "Nano Banana Pro not accessible"

**Solution**: Verify scientific-slides skill is installed

```bash
# Check skill directory
ls skills/scientific-slides/

# Test script manually
python skills/scientific-slides/scripts/generate_slide_image_ai.py --help
```

---

## CLI Reference

```bash
scrna-agent report --help
```

**Options:**
- `--output-dir, -o` (required): Directory with analysis results
- `--title, -t`: Presentation title (default: "scRNA-seq Analysis Results")

**Examples:**

```bash
# Basic
scrna-agent report -o ./outputs

# With custom title
scrna-agent report -o ./outputs -t "Multi-Tissue Atlas Study"

# Full workflow
scrna-agent annotate --file data.h5ad --tissue "Spleen"
scrna-agent report -o ./celltypist_outputs
```

---

## Output Files

After generation, you'll find:

```
celltypist_outputs/
├── celltypist_predictions.csv      # Analysis results
├── tissue_model_usage.csv          # Models used
├── run_metadata.json               # Metadata
├── presentation.pdf                # ← Final presentation
└── slides/
    ├── slide_01.png                # Individual slides
    ├── slide_02.png
    ├── ...
    └── slide_08.png
```

---

## Best Practices

1. **Run annotation first**: Always generate analysis results before creating presentation

2. **Use descriptive titles**: Make titles specific to your study
   ```bash
   scrna-agent report -o ./outputs -t "Multi-Tissue B Cell Atlas"
   ```

3. **Check results**: Verify analysis results are complete before generating

4. **Review slides**: Check individual slides in `slides/` directory

5. **Iterate**: Regenerate with different titles/customizations as needed

---

## Advanced Usage

### Batch Processing

Generate presentations for multiple datasets:

```bash
#!/bin/bash
for dir in outputs/sample_*; do
    sample=$(basename $dir)
    scrna-agent report \
      -o "$dir" \
      -t "Analysis: $sample"
done
```

### Custom Templates

Create your own slide generator:

```python
from scrna_agent.tools.report_generator import PresentationGenerator

class CustomGenerator(PresentationGenerator):
    def create_slide_plan(self, results, title):
        # Add your custom slides
        slides = super().create_slide_plan(results, title)

        # Add extra slide
        slides.append({
            "type": "content",
            "title": "My Custom Analysis",
            "prompt": "Create slide showing..."
        })

        return slides

# Use custom generator
gen = CustomGenerator("./outputs")
pdf = gen.generate_presentation("Custom Study")
```

---

## FAQ

**Q: Can I edit the PDF after generation?**
A: Yes, use PDF editors like Adobe Acrobat or online tools. Or regenerate with custom slides.

**Q: Can I change the design/colors?**
A: Yes, edit the prompts in `report_generator.py` to specify different colors/styles.

**Q: How long does generation take?**
A: ~1-2 minutes per slide. 8 slides = ~10-15 minutes total.

**Q: Can I add my own slides?**
A: Yes, create custom slide specifications and add to the plan.

**Q: Does this require internet?**
A: Yes, Nano Banana Pro AI requires internet connection.

---

## Examples

See `examples/generate_presentation.py` for:
- Simple generation
- Custom control
- Complete workflows

---

## Support

For issues:
1. Check this guide
2. Verify input files exist
3. Test Nano Banana Pro manually
4. Check CLI help: `scrna-agent report --help`

---

## Summary

**Installation**: Skill already included in scrna_agent/skills/

**Usage**:
```bash
scrna-agent report -o ./outputs -t "My Study"
```

**Output**: Professional PDF presentation ready for presenting

**Time**: ~10-15 minutes for 8-slide deck

**Customization**: Edit prompts in report_generator.py or create custom slides

---

**Ready to present your science!** 🎉
