"""
Report Generator for scRNA-seq Analysis
Generates PPT presentations from analysis results
"""

import os
import json
import subprocess
from pathlib import Path
from typing import Dict, List, Any, Optional
import pandas as pd


class PresentationGenerator:
    """Generate PowerPoint presentations from scRNA-seq analysis results."""

    def __init__(self, output_dir: str):
        """
        Initialize presentation generator.

        Args:
            output_dir: Directory containing analysis results
        """
        self.output_dir = Path(output_dir)
        self.slides_dir = self.output_dir / "slides"
        self.slides_dir.mkdir(exist_ok=True)

        # Get skills directory
        package_dir = Path(__file__).parent.parent
        self.skills_dir = package_dir / "skills" / "scientific-slides"
        self.scripts_dir = self.skills_dir / "scripts"

    def load_analysis_results(self) -> Dict[str, Any]:
        """Load analysis results from output directory."""
        results = {}

        # Load tissue model usage
        tissue_usage_file = self.output_dir / "tissue_model_usage.csv"
        if tissue_usage_file.exists():
            results["tissue_usage"] = pd.read_csv(tissue_usage_file)

        # Load predictions
        predictions_file = self.output_dir / "celltypist_predictions.csv"
        if predictions_file.exists():
            results["predictions"] = pd.read_csv(predictions_file)

        # Load metadata
        metadata_file = self.output_dir / "run_metadata.json"
        if metadata_file.exists():
            with open(metadata_file) as f:
                results["metadata"] = json.load(f)

        return results

    def create_slide_plan(self, results: Dict[str, Any], title: str = "scRNA-seq Analysis Results") -> List[Dict[str, str]]:
        """
        Create slide plan from analysis results.

        Args:
            results: Analysis results dictionary
            title: Presentation title

        Returns:
            List of slide specifications
        """
        slides = []

        # Slide 1: Title
        slides.append({
            "type": "title",
            "title": title,
            "subtitle": "Tissue-Aware Cell Type Annotation",
            "author": "scRNA-agent",
            "prompt": f"Create a professional title slide with title '{title}', subtitle 'Tissue-Aware Cell Type Annotation', author 'scRNA-agent'. Use a modern scientific design with dark blue background, white text, and subtle molecular/cell imagery in the background."
        })

        # Slide 2: Overview
        if "metadata" in results:
            metadata = results["metadata"]
            total_cells = metadata.get("total_cells", "N/A")
            total_tissues = metadata.get("total_tissues", "N/A")

            slides.append({
                "type": "content",
                "title": "Analysis Overview",
                "prompt": f"Create a slide titled 'Analysis Overview' with three key statistics displayed prominently: '\\n- Total Cells: {total_cells:,}\\n- Tissues Analyzed: {total_tissues}\\n- Annotation Method: CellTypist (Tissue-Aware)\\n\\nUse large numbers, icons for each stat, clean modern design with dark blue background."
            })

        # Slide 3: Tissue Distribution
        if "tissue_usage" in results:
            tissue_df = results["tissue_usage"]
            tissue_info = []
            for _, row in tissue_df.iterrows():
                tissue_info.append(f"- {row['tissue']}: {row['n_cells']:,} cells")

            tissue_text = "\\n".join(tissue_info[:5])  # Top 5 tissues

            slides.append({
                "type": "content",
                "title": "Tissue Distribution",
                "prompt": f"Create a slide titled 'Tissue Distribution' showing:\\n{tissue_text}\\n\\nUse a horizontal bar chart or visual representation. Modern design, dark blue background, white text with gold accents for numbers."
            })

        # Slide 4: Model Selection Strategy
        slides.append({
            "type": "content",
            "title": "Model Selection Strategy",
            "prompt": "Create a slide titled 'Model Selection Strategy' showing the priority hierarchy:\\n\\n1. Override (tissue|sample_type) → Exact match\\n2. Tissue → Tissue-specific model\\n3. Sample Type → Sample type-specific\\n4. Default → Fallback model\\n\\nUse a flowchart or decision tree visualization. Dark blue background, clean modern design."
        })

        # Slide 5: Models Used
        if "tissue_usage" in results:
            tissue_df = results["tissue_usage"]
            model_info = []
            for _, row in tissue_df.iterrows():
                model_info.append(f"- {row['tissue']}: {row['model_used']} ({row['selection_rule']})")

            model_text = "\\n".join(model_info[:5])

            slides.append({
                "type": "content",
                "title": "Models Used per Tissue",
                "prompt": f"Create a slide titled 'Models Used per Tissue' showing:\\n{model_text}\\n\\nUse a table or list format with icons. Dark blue background, modern scientific design."
            })

        # Slide 6: Cell Type Distribution
        if "predictions" in results:
            pred_df = results["predictions"]
            top_celltypes = pred_df["label"].value_counts().head(10)
            celltype_info = []
            for celltype, count in top_celltypes.items():
                celltype_info.append(f"- {celltype}: {count:,} cells")

            celltype_text = "\\n".join(celltype_info)

            slides.append({
                "type": "content",
                "title": "Top Cell Types Identified",
                "prompt": f"Create a slide titled 'Top Cell Types Identified' showing:\\n{celltype_text}\\n\\nUse a pie chart or donut chart visualization. Dark blue background, colorful segments, modern design."
            })

        # Slide 7: Confidence Scores
        if "predictions" in results:
            pred_df = results["predictions"]
            mean_conf = pred_df["confidence"].mean()
            median_conf = pred_df["confidence"].median()
            low_conf = (pred_df["confidence"] < 0.5).sum()

            slides.append({
                "type": "content",
                "title": "Prediction Confidence",
                "prompt": f"Create a slide titled 'Prediction Confidence' showing:\\n- Mean Confidence: {mean_conf:.3f}\\n- Median Confidence: {median_conf:.3f}\\n- Low Confidence (<0.5): {low_conf:,} cells\\n\\nUse gauge charts or progress bars. Dark blue background, gold/green for good confidence, red for low."
            })

        # Slide 8: Summary
        slides.append({
            "type": "content",
            "title": "Summary",
            "prompt": "Create a slide titled 'Summary' with key takeaways:\\n- Successfully annotated cells across multiple tissues\\n- Automatic model selection based on tissue type\\n- High confidence predictions\\n- Ready for downstream analysis\\n\\nUse checkmarks and clean bullet points. Dark blue background, professional design."
        })

        return slides

    def generate_slide(self, slide_spec: Dict[str, str], slide_num: int, previous_slide: Optional[Path] = None) -> Path:
        """
        Generate a single slide using Nano Banana Pro.

        Args:
            slide_spec: Slide specification with prompt
            slide_num: Slide number
            previous_slide: Path to previous slide for style consistency

        Returns:
            Path to generated slide image
        """
        output_file = self.slides_dir / f"slide_{slide_num:02d}.png"

        # Build command
        cmd = [
            "python",
            str(self.scripts_dir / "generate_slide_image_ai.py"),
            slide_spec["prompt"],  # Prompt is a positional argument
            "-o", str(output_file)
        ]

        # Attach previous slide for style consistency
        if previous_slide and previous_slide.exists():
            cmd.extend(["--attach", str(previous_slide)])

        print(f"Generating slide {slide_num}: {slide_spec.get('title', 'Untitled')}...")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if result.returncode != 0:
                print(f"Warning: Slide generation had issues: {result.stderr}")
        except subprocess.TimeoutExpired:
            print(f"Warning: Slide {slide_num} generation timed out")
        except Exception as e:
            print(f"Warning: Error generating slide {slide_num}: {e}")

        return output_file

    def combine_slides_to_pdf(self, output_pdf: Path) -> None:
        """
        Combine all slides into a PDF.

        Args:
            output_pdf: Path to output PDF file
        """
        # Get all slide images
        slide_images = sorted(self.slides_dir.glob("slide_*.png"))

        if not slide_images:
            raise ValueError("No slides found to combine")

        cmd = [
            "python",
            str(self.scripts_dir / "slides_to_pdf.py"),
            str(output_pdf),
            *[str(img) for img in slide_images]
        ]

        print(f"Combining {len(slide_images)} slides into PDF...")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            if result.returncode != 0:
                print(f"Warning: PDF generation had issues: {result.stderr}")
            else:
                print(f"✓ PDF created: {output_pdf}")
        except Exception as e:
            print(f"Error creating PDF: {e}")

    def generate_presentation(self, title: str = "scRNA-seq Analysis Results") -> Path:
        """
        Generate complete presentation from analysis results.

        Args:
            title: Presentation title

        Returns:
            Path to generated PDF
        """
        print("=" * 60)
        print("Generating Presentation")
        print("=" * 60)

        # Load results
        print("Loading analysis results...")
        results = self.load_analysis_results()

        if not results:
            raise ValueError("No analysis results found in output directory")

        # Create slide plan
        print("Creating slide plan...")
        slides = self.create_slide_plan(results, title)
        print(f"Planned {len(slides)} slides")

        # Generate each slide
        previous_slide = None
        for i, slide_spec in enumerate(slides, 1):
            slide_path = self.generate_slide(slide_spec, i, previous_slide)
            previous_slide = slide_path

        # Combine to PDF
        output_pdf = self.output_dir / "presentation.pdf"
        self.combine_slides_to_pdf(output_pdf)

        print("=" * 60)
        print(f"✓ Presentation complete: {output_pdf}")
        print("=" * 60)

        return output_pdf


def generate_report(
    output_dir: str,
    title: str = "scRNA-seq Analysis Results"
) -> str:
    """
    Generate PPT presentation from analysis results.

    Args:
        output_dir: Directory containing analysis results
        title: Presentation title

    Returns:
        Path to generated PDF
    """
    generator = PresentationGenerator(output_dir)
    pdf_path = generator.generate_presentation(title)
    return str(pdf_path)
