# Neuronal-silence-as-an-epileptic-seizure-biomarker

This repository contains the codes and resources associated with the manuscript *"Neuronal silence as an epileptic seizure biomarker"*.  

In this work, we propose and test a **new potential biomarker** for seizure prediction and detection in neuronal model time series. Our approach focuses on the **silence duration** of a neuronal population, showing that **prolonged silent periods often precede abnormal synchronization events**, a feature of seizure onset.  

The repository is divided into two main components:  
1. **Neuronal Model Simulation (C code)** – reproduces simulation results from the paper, specifically *Figure 5(a)*.  
2. **Random Forest Algorithm (Python)** – implements our machine learning analysis using a pre-trained Random Forest model and an example silence time series.  

The instructions on how to run the code and modify parameters are described below.

> **Citation:** If you use this repository in your work, please cite our paper (DOI to be added upon publication).  

---

## Neuronal Model Simulation

All simulations in the paper were implemented in **C**. In this repository, we provide one example code that generates *Figure 5(a)* from the manuscript.  

**Requirements:**  
- A C compiler such as `gcc` or `icc`  

**Running the code:**  
Executing the simulation produces three `.dat` files:  
- **Kuramoto order parameter** `R(t)` – measures global synchronization  
- **Coefficient of variation** `CV(t)` – quantifies spike irregularity  
- **Mean silence time** `<T>` – the biomarker proposed in this work  

By default, the code simulates **250 seconds** of network activity. You can adjust:  
- **Simulation time** (shorter or longer runs)  
- **Coupling strength**  
- **Other model parameters**  

This flexibility allows exploration of how `<T>` behaves under different dynamical regimes.  

**Visualization:**  
A Jupyter Notebook (`plot.ipynb`) is included for plotting and analyzing the results.  

---

## Predictive Random Forest Algorithm

This folder contains the Python notebook implementation of the **predictive Random Forest** model used in our study to identify seizure-related patterns from the mean silence time series.  

In our study, the Random Forest was trained to **predict the Kuramoto Order Parameter** `R(t)` from the `<T>` (mean silence duration) time series.  
To achieve this, we provide the model with **two time-delayed values** of `<T>`, corresponding to a total delay of **10 ms**, so that the algorithm must **predict the next point** of `R(t)`.  
This design forces the model to capture short-term temporal dependencies and use them for forward prediction.  

This predictive model demonstrated strong performance in estimating the **synchronization level** of the neuronal network across different stimulation frequencies and coupling strengths, while using the **same training set**. This highlights the robustness and generalization capacity of our approach.  

**Contents of this folder:**  
- Four **example mean silence time series** from the simulations  
- Python code to:
  - Train the model 
  - Perform predictions on new data  
  - Evaluate performance using standard classification metrics  

**What you can do with it:**  
- **Reproduce the predictive results** shown in the paper  
- **Test the model with your own time series** of mean silence duration  
- **Retrain** the Random Forest using new datasets or different parameters  

By running this code, you can see how our biomarker `<T>` can be used in real time to provide the synchronization level.

