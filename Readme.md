# 🧬 Drug-Protein Interaction Predictor (BioPredictor V3.6)

![Python](https://img.shields.io/badge/Python-3.9%2B-blue)
![Flask](https://img.shields.io/badge/Framework-Flask-green)
![ML](https://img.shields.io/badge/AI-Random_Forest-orange)
![Bioinformatics](https://img.shields.io/badge/Bio-RDKit-red)

> **An AI-powered bioinformatics tool that accelerates early-stage drug discovery by predicting binding affinity and visualizing molecular interactions in 3D.**

---

## 📌 Overview
Developing a new drug takes 10–15 years and costs billions. A major bottleneck is screening millions of compounds to find one that binds to a disease target.

**BioPredictor V3.6** solves this by using Machine Learning to instantly predict the interaction probability between a drug and a protein. It combines a **Random Forest Backend** with a **WebGL Frontend** to provide both mathematical confidence and visual proof of the "Lock-and-Key" mechanism.

---

## 🚀 Key Features

### 1. 🧠 The AI Prediction Engine
- **Dual-Input Processing:** Fuses chemical data (Drugs) and biological data (Proteins) into a single feature vector.
- **Algorithms:** Uses **Morgan Fingerprints** for drug substructure analysis and **Amino Acid Composition (AAC)** for protein profiling.
- **Scientific Reality Filter:** A specialized heuristic that normalizes confidence scores to mimic biological uncertainty (preventing unrealistic 100% predictions).

### 2. 🔬 Interactive 3D Visualization
- **Real-Time Rendering:** Uses `3Dmol.js` to render hardware-accelerated protein structures in the browser.
- **Database Integration:** Fetches verified crystallographic data directly from the **RCSB Protein Data Bank (PDB)**.
- **Visual Modes:** Supports Ribbon, Stick, and Sphere representations to analyze binding pockets.

---

## 🛠️ Technical Architecture

| Component | Technology Used | Role |
| :--- | :--- | :--- |
| **Backend** | Python, Flask | REST API & Server Logic |
| **AI Model** | Scikit-Learn (Random Forest) | Binding Probability Inference |
| **Cheminformatics** | RDKit | SMILES parsing & Morgan Fingerprints |
| **Data Processing** | Pandas, NumPy | Vectorization & Matrix Math |
| **Frontend** | HTML5, CSS3, JavaScript | UI & Interaction Logic |
| **Visualization** | 3Dmol.js (WebGL) | Rendering PDB files |

---

## ⚙️ Installation & Setup

1. **Clone the repository**
   ```bash
   git clone [https://github.com/huzaifashah78/DRUG-PROTIEN-INTERACTION-PREDICTOR.git](https://github.com/huzaifashah78/DRUG-PROTIEN-INTERACTION-PREDICTOR.git)
   cd DRUG-PROTIEN-INTERACTION-PREDICTOR