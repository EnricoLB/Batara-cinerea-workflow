
## Platform-Agnostic SAM3 + BioEncoder Workflow

### Requirements
- Python 3.10
- A machine with a GPU (or Google Colab with GPU runtime)
- I used HiperGator. The workflow for HiperGator is in the same folder as this file on GitHub.

---

### Step 1 — Create and activate environment

```bash
# Using conda (local machine or any HPC)
conda create -n sam3_env python=3.10 -y
conda activate sam3_env

# OR using venv (no conda needed)
python3.10 -m venv sam3_env
source sam3_env/bin/activate        # Mac/Linux
sam3_env\Scripts\activate           # Windows
```

---

### Step 2 — Install dependencies

```bash
pip install autodistill autodistill-sam3 roboflow inference \
            opencv-python numpy scikit-learn sam3 bioencoder
```

---

### Step 3 — Set Roboflow API key

```bash
# Mac/Linux
export ROBOFLOW_API_KEY="YOUR_KEY"

# Windows
set ROBOFLOW_API_KEY="YOUR_KEY"
```

To get your key: go to https://app.roboflow.com/settings/api

---

### Step 4 — Organize files (optional)

```bash
# Create dataset structure BioEncoder expects
mkdir -p dataset/cinerea dataset/argentina dataset/excubitor

# Copy segmented images into class folders
cp segmented_cinerea/*   dataset/cinerea/
cp segmented_argentina/* dataset/argentina/
cp segmented_excubitor/* dataset/excubitor/
```

---

### Step 5 — Run SAM3 segmentation

```python
# sam3_segmentation.py
from autodistill_sam3 import SAM3
from autodistill.detection import CaptionOntology
import os

model = SAM3()

input_dir  = "images/"       # folder with your raw images
output_dir = "segmented/"    # where masked images will be saved
os.makedirs(output_dir, exist_ok=True)

for img_file in os.listdir(input_dir):
    if img_file.endswith((".jpg", ".png", ".jpeg")):
        img_path = os.path.join(input_dir, img_file)
        results  = model.predict(img_path)
        results.save(os.path.join(output_dir, img_file))
```

```bash
python sam3_segmentation.py
```

---

### Step 6 — Configure BioEncoder workspace

```python
import bioencoder
bioencoder.configure(root_dir="bioencoder_wd", run_name="run_v1")
```

---

### Step 7 — Split dataset

```bash
bioencoder_split_dataset --image-dir "dataset" --max-ratio 6 --random-seed 42
```

---

### Step 8 — ⚠️ Fix number of classes (critical step)

Edit both YAML files before training [3]:

```bash
nano yml_files/train_stage2.yml
nano yml_files/swa_stage2.yml
```

Set `num_classes` and `classes` to match your number of folders (e.g. 3 for *cinerea*, *argentina*, *excubitor*).

---

### Step 9 — Train model

```bash
bioencoder_train --config-path "yml_files/train_stage1.yml"
bioencoder_swa   --config-path "yml_files/swa_stage1.yml"
bioencoder_interactive_plots --config-path "yml_files/plot_stage1.yml"

bioencoder_train --config-path "yml_files/train_stage2.yml" --overwrite
bioencoder_swa   --config-path "yml_files/swa_stage2.yml"   --overwrite
```

---

### Step 10 — GradCAM model explorer

```bash
# Fix known img_size config bug first
echo "img_size: 224" > temp.yml
cat yml_files/explore_stage2.yml >> temp.yml
mv temp.yml yml_files/explore_stage2.yml

sed -i 's/img_size: \[224, 224\]/img_size: 224/g' yml_files/explore_stage2.yml

# Launch
bioencoder_model_explorer --config-path "yml_files/explore_stage2.yml"
```

Then open `http://localhost:8501` in your browser.

---

### Google Colab alternative

If you have no local GPU, the entire workflow above runs on Colab. Just add this at the top of your notebook:

```python
from google.colab import drive
drive.mount("/content/drive")

# Then replace all local paths like "dataset/" with
# "/content/drive/MyDrive/your_folder/dataset/"
```

And replace the `conda`/`venv` setup with:

```python
!pip install autodistill autodistill-sam3 roboflow inference \
             opencv-python numpy scikit-learn sam3 bioencoder
```
