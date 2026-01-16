# PPT Generation - Quick Start

## 빠른 시작 (한국어)

### 1단계: 분석 실행

```bash
# Tissue 매핑 파일 생성
cat > tissue_map.csv << EOF
sample,tissue
sample1,Spleen
sample2,Blood
sample3,Liver
EOF

# Annotation 실행
scrna-agent annotate \
  --file your_data.h5ad \
  --tissue-map tissue_map.csv \
  --output annotated.h5ad
```

### 2단계: PPT 생성

```bash
# PPT 자동 생성
scrna-agent report \
  --output-dir ./celltypist_outputs \
  --title "내 scRNA-seq 연구"
```

✅ **완료!** `./celltypist_outputs/presentation.pdf` 파일이 생성됩니다.

---

## Quick Start (English)

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
  --file your_data.h5ad \
  --tissue-map tissue_map.csv \
  --output annotated.h5ad
```

### Step 2: Generate PPT

```bash
# Auto-generate presentation
scrna-agent report \
  --output-dir ./celltypist_outputs \
  --title "My scRNA-seq Study"
```

✅ **Done!** Find your presentation at `./celltypist_outputs/presentation.pdf`

---

## Python API

```python
# Run annotation
from scrna_agent.tools.annotation_tools import *

load_adata_for_annotation("data.h5ad")
add_tissue_metadata(tissue_map="map.csv")
celltypist_annotate_tissue_aware(output_dir="./outputs")

# Generate PPT
from scrna_agent.tools.report_generator import generate_report

pdf = generate_report("./outputs", "My Study")
print(f"Presentation: {pdf}")
```

---

## What's Included

자동 생성되는 슬라이드:

1. **Title Slide** - 제목, 저자
2. **Overview** - 전체 통계
3. **Tissue Distribution** - 조직별 세포 수
4. **Model Selection** - 모델 선택 전략
5. **Models Used** - 사용된 모델
6. **Cell Types** - 발견된 세포 타입
7. **Confidence** - 예측 신뢰도
8. **Summary** - 요약

---

## Requirements

```bash
# Install dependencies
pip install pillow reportlab

# Verify skill is installed
ls scrna_agent/skills/scientific-slides/
```

---

## Troubleshooting

**문제**: "No results found"
**해결**: Annotation을 먼저 실행하세요

```bash
scrna-agent annotate --file data.h5ad --tissue "Spleen"
scrna-agent report -o ./celltypist_outputs
```

**문제**: "Nano Banana Pro error"
**해결**: 인터넷 연결 확인 및 스킬 디렉토리 확인

```bash
ls scrna_agent/skills/scientific-slides/
```

---

## Full Documentation

자세한 내용은:
- **Full Guide**: `docs/PPT_GENERATION_GUIDE.md`
- **Examples**: `examples/generate_presentation.py`

---

## Quick Commands

```bash
# Help
scrna-agent report --help

# Basic usage
scrna-agent report -o ./outputs

# With custom title
scrna-agent report -o ./outputs -t "Multi-Tissue Atlas"

# Complete workflow
scrna-agent annotate --file data.h5ad --tissue "Spleen" && \
scrna-agent report -o ./celltypist_outputs
```

---

**시작하세요!** 🚀 Just run annotation, then `scrna-agent report`!
