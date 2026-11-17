# Micro-PIV Blood Flow Analysis (Hemodynamics & Biomechanics)

A quantitative Particle Image Velocimetry (PIV) analysis pipeline for measuring in-plane velocity fields, strain rate tensors, and viscous stress distributions in blood micro-flow.

---

## 📌 Project Overview
This project processes consecutive micro-PIV image pairs to quantify hemodynamic parameters in microcirculation.

- **Image Pair:** `image1.bmp` and `Image2.bmp` ($1280 \times 1024$ pixels)
- **Time Separation ($\Delta t$):** $300\text{ ms} = 0.300\text{ s}$
- **Spatial Resolution:** $1\text{ pixel} = 1.4\ \mu\text{m}$ ($1.4 \times 10^{-6}\text{ m}$) in both $x$ and $y$
- **Working Fluid:** Whole human blood ($\mu = 0.0035\text{ Pa}\cdot\text{s}$, $\rho = 1060\text{ kg/m}^3$)

---

## 🔬 Mathematical Formulations

### 1. Velocity Field ($U, V, W$)
Displacements $dx$ and $dy$ are determined via FFT-based cross-correlation with 3-point Gaussian sub-pixel peak fitting:
$$U = \frac{dx \cdot s}{\Delta t}, \quad V = -\frac{dy \cdot s}{\Delta t}, \quad W = 0$$
$$\text{Speed} = \sqrt{U^2 + V^2 + W^2}$$
*(Where $s = 1.4\ \mu\text{m/pixel}$, and $W = 0$ due to planar 2D PIV).*

### 2. Rate-of-Deformation (Strain-Rate) Tensor ($E_{ij}$)
Calculated using second-order central spatial differences:
$$E_{xx} = \frac{\partial U}{\partial x}, \quad E_{yy} = \frac{\partial V}{\partial y}$$
$$E_{xy} = \frac{1}{2}\left(\frac{\partial U}{\partial y} + \frac{\partial V}{\partial x}\right)$$
$$\|E\| = \sqrt{E_{xx}^2 + E_{yy}^2 + 2E_{xy}^2}$$

### 3. Viscous Stress Tensor ($\tau_{ij}$)
Modeled using the incompressible Newtonian constitutive relationship:
$$\tau_{ij} = 2\mu E_{ij}$$
$$\tau_{xx} = 2\mu E_{xx}, \quad \tau_{yy} = 2\mu E_{yy}$$
$$\tau_{xy} = 2\mu E_{xy} = \mu\left(\frac{\partial U}{\partial y} + \frac{\partial V}{\partial x}\right)$$
$$\|\tau\| = \sqrt{\tau_{xx}^2 + \tau_{yy}^2 + 2\tau_{xy}^2}$$

---

## 📊 Summary of Results

| Parameter | Median | Mean | Peak Value | Unit |
| :--- | :--- | :--- | :--- | :--- |
| **Horizontal Velocity ($U$)** | -0.030 | -0.033 | -0.062 | mm/s |
| **Vertical Velocity ($V$)** | 0.006 | 0.006 | 0.055 | mm/s |
| **Flow Speed ($|V|$)** | 0.036 | 0.038 | 0.063 | mm/s |
| **Shear Strain Rate ($E_{xy}$)** | -0.030 | -0.032 | -0.069 | s⁻¹ |
| **Strain Rate Magnitude** | 0.052 | 0.053 | 0.135 | s⁻¹ |
| **Viscous Shear Stress ($\tau_{xy}$)** | -0.213 | -0.224 | -0.484 | mPa |
| **Viscous Stress Magnitude** | 0.364 | 0.373 | 0.945 | mPa |

---

## 📁 Repository Structure
```
├── image1.bmp                 # Frame 1 at time t
├── Image2.bmp                 # Frame 2 at time t + 300 ms
├── PIV_Blood_Flow_Analysis.py # Main PIV analysis script
├── README.md                  # Project documentation
└── PIV_Blood_Flow_Output/     # Output results directory
    ├── PIV_Blood_Flow_Results.csv
    ├── velocity_vectors.png
    ├── velocity_magnitude.png
    ├── U_velocity.png
    ├── V_velocity.png
    ├── W_assumed_zero.png
    ├── strain_magnitude.png
    ├── Exx.png, Eyy.png, Exy.png
    ├── stress_magnitude.png
    └── tau_xx.png, tau_yy.png, tau_xy.png
```

---

## 🚀 How to Run
```bash
# Clone the repository
git clone https://github.com/santhosh1728/biomech_project_2_PIV-Blood-Flow-Analysis.git
cd biomech_project_2_PIV-Blood-Flow-Analysis

# Install dependencies
pip install numpy pandas scipy matplotlib pillow opencv-python

# Run analysis
python PIV_Blood_Flow_Analysis.py
```
Outputs are automatically saved into `PIV_Blood_Flow_Output/`.
