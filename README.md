# Multiscale Modeling of Materials Using Machine Learning

## Project Overview
This project focuses on **optimizing the microstructural properties of copper heat pipes** to enhance their performance. By integrating **machine learning (ML)-derived constitutive models** with finite element analysis (FEA), we explored the impact of grain size on stress handling and thermal efficiency. The concept of **multiscale modeling** was central to connecting the grain-scale microstructure to the larger, component-level material behavior.

---

## Objectives
1. **Analyze the effect of grain size** on the strength and heat resistance of copper heat pipes.
2. Develop a **custom user material (UMAT)** subroutine based on an ML-derived stress-strain relationship.
3. **Validate** the UMAT model against experimental and simulation results.

---

## Methodology

### 1. Multiscale Modeling
- **Micro to Macro Connection**:
  - Simulated **grain-scale behavior** using the VPSC (Visco-Plastic Self-Consistent) code.
  - Integrated grain-scale data (stress, strain rate, grain size) into the **macroscale FEA** model using the UMAT subroutine.
  
- **Purpose**:
  - To understand how microstructural features like grain size influence macroscale material properties (e.g., stress, strain, and thermal conductivity).

### 2. Data Generation
- **VPSC Simulation Setup**: 
  - Generated 76 datasets using `vpsc7.exe` for **19 grain sizes** (5–95 µm) at **4 strain rates** (10⁻² to 10⁻⁵).
  - Each dataset contained stress-strain data for different combinations of grain size and strain rate.

- **Dataset Concatenation**: Combined all datasets for ML training and analysis.

### 3. Machine Learning Model
- **Polynomial Regression (Degree 2)**:
  - Modeled **stress** as a function of **strain**, **strain rate**, and **grain size**.
  - Derived the following stress equation:
    ```text
    Stress = 115.87 + 581.17 * Strain + 24694.13 * Strain Rate + 1168663.06 * Grain Size 
             - 1368.07 * Strain² + 12543.13 * Strain * Strain Rate 
             - 29095.13 * Strain * Grain Size - 2102507.11 * Strain Rate² 
             - 9820971.89 * Strain Rate * Grain Size + 8222685370.00 * Grain Size²
    ```
  - Achieved **R² = 0.948**.

### 4. UMAT Code Development
- Integrated the **ML-derived stress equation** into the UMAT subroutine in **Abaqus**.
- Differentiated the equation to compute:
  - **Stress Flow (sf)**
  - **Strain Tensor Update**

- Enhanced the UMAT to include **plastic strain adjustments** for dynamic stress-strain calculations.

### 5. Validation
- Compared **stress-strain plots** from the UMAT model against VPSC simulation results.
- Result: **Close match**, validating the ML-based approach.

---

## Simulation Details
- **Component**: Hollow copper pipe.
  - Inner Diameter: 8 mm
  - Outer Diameter: 10 mm
  - Length: 100 mm

- **FEA Setup**:
  - Applied **radial load** to simulate real-world conditions.
  - Two scenarios:
    1. **With Microstructure**: Grain size incorporated using the UMAT.
    2. **Without Microstructure**: Uniform material response.

---

## Multiscale Insights
- **Grain Size Effect**:
  - Smaller grain sizes increased stress due to the **Hall-Petch effect**, improving material strength.
  - Enhanced the **strain rate sensitivity**, reflecting realistic material behavior under varying loads.

- **Integration of Scales**:
  - By incorporating grain-scale behavior into the macroscale FEA model, the simulation accurately captured localized stress variations.
  - The multiscale approach highlighted the **interdependence of microstructure and component performance**, bridging the gap between material science and engineering applications.

---

## Key Findings
1. **Smaller grain sizes** increased material strength due to the Hall-Petch effect.
2. Incorporating microstructure led to **localized stress variations**, improving accuracy in simulations.
3. **Thermal Performance**:
   - Smaller grains reduced heat resistance, enabling faster heat transfer without damaging the pipe.

---

## Conclusion
The integration of microstructural details like grain size into FEA simulations:
- Improved stress and heat transfer predictions.
- Highlighted the importance of **multiscale modeling** in optimizing material performance.

This study demonstrates the potential of combining machine learning with FEM to achieve realistic and efficient material design.

---

## Technologies and Tools Used
- **Programming**: Fortran
- **Simulation**: Abaqus, VPSC
- **Machine Learning**: Python (Polynomial Regression)
- **Software**: Ansys for custom material definition

---

## Visual Results
### Stress-Strain Curves
1. **With Microstructure**: 
   - Stress levels significantly increased with smaller grain sizes.

2. **Without Microstructure**: 
   - Uniform stress distribution observed.

---

Feel free to explore the code and dataset in this repository for further insights.
