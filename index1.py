import os
import pandas as pd
import numpy as np
from scipy.stats import linregress
from scipy.interpolate import make_interp_spline
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns

# Ensure the /img directory exists
folder = "img1/"
os.makedirs(folder, exist_ok=True)

# Set Seaborn style for better aesthetics
sns.set(style="whitegrid")

# Calibration Curve Data (Potassium Dichromate - K2Cr2O7)
calibration_data = {
    "Absorbance": [0, 0.15, 0.35, 0.69, 0.95, 0.98],  # Absorbance Values
    "Concentration (mg/mL)": [0, 0.2, 0.4, 0.6, 0.8, 1]  # Concentration Values
}

# Convert to DataFrame
calib_df = pd.DataFrame(calibration_data)

# Perform Linear Regression (y = mx + b)
slope, intercept, r_value, _, _ = linregress(calib_df["Concentration (mg/mL)"], calib_df["Absorbance"])

# Print Calibration Equation
print(f"Calibration Curve Equation: Absorbance = ({slope:.4f} * Concentration) + {intercept:.4f}")
print(f"R-squared: {r_value**2:.4f}")

# Define x and y values
x = np.array(calib_df["Concentration (mg/mL)"])
y = np.array(calib_df["Absorbance"])

# Generate Smooth Curve (Cubic Spline)
x_smooth = np.linspace(x.min(), x.max(), 200)  # More points for smooth curve
spline = make_interp_spline(x, y, k=3)  # k=3 for cubic spline
y_smooth = spline(x_smooth)

# --- Plot 1: Scatter + Linear Fit ---
plt.figure(figsize=(8, 5))
sns.scatterplot(x=x, y=y, color='blue', label="Calibration Data")
plt.plot(x, slope * x + intercept, color='red', linestyle='--', label=f"Fit: y={slope:.4f}x + {intercept:.4f}")
plt.text(0.1, 0.8, f"R² = {r_value**2:.4f}", fontsize=12, color="black", transform=plt.gca().transAxes)
plt.xlabel("Concentration (mg/mL)")
plt.ylabel("Absorbance")
plt.title("Calibration Curve: Absorbance vs. Concentration")
plt.legend()
plt.savefig(f"{folder}/calibration_linear (K2Cr2O7).png", dpi=300, bbox_inches="tight")
# plt.show()

# --- Plot 2: Scatter + Smooth Curve ---
plt.figure(figsize=(8, 5))
sns.scatterplot(x=x, y=y, color='blue', label="Calibration Data")
plt.plot(x_smooth, y_smooth, color='black', linestyle='-', linewidth=2, label="Cubic Spline Fit")
plt.plot(x, slope * x + intercept, color='red', linestyle='--', label=f"Fit: y={slope:.4f}x + {intercept:.4f}")
plt.text(0.1, 0.8, f"R² = {r_value**2:.4f}", fontsize=12, color="black", transform=plt.gca().transAxes)
plt.xlabel("Concentration (mg/mL)")
plt.ylabel("Absorbance")
plt.title("Calibration Curve (K2Cr2O7)")
plt.legend()
plt.savefig(f"{folder}/calibration_smooth (K2Cr2O7).png", dpi=300, bbox_inches="tight")
# plt.show()

# Cell growth 
# Experimental Data (time, absorbance, pH, temp)
experimental_data = {
    'Time (hrs)': [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0],  # Time values
    'Absorbance': [0.29, 0.32, 0.47, 0.60, 0.70, 0.74, 0.76, 0.77, 0.77, 0.70, 0.57, 0.54, 0.42]  # Absorbance values
}

cell_growth_df = pd.DataFrame(experimental_data)

# Convert Absorbance to Cell Growth Concentration using Calibration curve 
cell_growth_df["Cell Growth Concentration (mg/mL)"] = (cell_growth_df["Absorbance"] - intercept) / slope

# Initialize the first value as NaN since we can't compute µ at time 0
mu_values = [0]

for i in range(1, len(cell_growth_df)):
    X0 = cell_growth_df["Cell Growth Concentration (mg/mL)"].iloc[i-1]  # Previous time concentration
    Xt = cell_growth_df["Cell Growth Concentration (mg/mL)"].iloc[i]  # Current time concentration
    t = cell_growth_df["Time (hrs)"].iloc[i] - cell_growth_df["Time (hrs)"].iloc[i-1]  # Time interval

    mu = (np.log(Xt / X0)) / t  # Specific growth rate calculation
    mu_values.append(mu)

# Add the computed growth rates to the DataFrame
cell_growth_df["Specific Growth Rate (µ)"] = mu_values

# Save DataFrame to CSV for easy identification
cell_growth_df.to_csv("cell_growth_data.csv", index=False)
print("Data saved successfully!")

# Display the updated DataFrame
cell_growth_df.head(2)

# --- Specific Growth Rate Plot ---
plt.figure(figsize=(8, 5))
x_smooth = np.linspace(min(cell_growth_df["Time (hrs)"]), max(cell_growth_df["Time (hrs)"]), 300)
y_smooth = make_interp_spline(cell_growth_df["Time (hrs)"], cell_growth_df["Specific Growth Rate (µ)"])(x_smooth)

sns.lineplot(x=x_smooth, y=y_smooth, color='b', label="Specific Growth Rate (µ)")
sns.scatterplot(x=cell_growth_df["Time (hrs)"], y=cell_growth_df["Specific Growth Rate (µ)"], color='b', label="Data Points")
plt.xlabel("Time (hrs)")
plt.ylabel("Specific Growth Rate (µ) (per hour)")
plt.title("Specific Growth Rate Over Time")
plt.legend()
plt.savefig(os.path.join(folder, "specific_growth_rate_curve.png"), dpi=300, bbox_inches="tight")
# plt.show()

# --- Cell Growth Concentration Curve ---
plt.figure(figsize=(8, 5))
x_smooth = np.linspace(min(cell_growth_df["Time (hrs)"]), max(cell_growth_df["Time (hrs)"]), 300)
y_smooth = make_interp_spline(cell_growth_df["Time (hrs)"], cell_growth_df["Cell Growth Concentration (mg/mL)"])(x_smooth)

sns.lineplot(x=x_smooth, y=y_smooth, color='b', label="Cell Growth")
sns.scatterplot(x=cell_growth_df["Time (hrs)"], y=cell_growth_df["Cell Growth Concentration (mg/mL)"], color='b', label="Data Points")
plt.xlabel("Time (hrs)")
plt.ylabel("Cell Growth Concentration (mg/mL)")
plt.title("Cell Growth Concentration Curve")
plt.legend()
plt.savefig(os.path.join(folder, "cell_growth_curve.png"), dpi=300, bbox_inches="tight")
# plt.show()

# Load experimental data
time = cell_growth_df["Time (hrs)"]
cell_concentration = cell_growth_df["Cell Growth Concentration (mg/mL)"]

# Smooth the curve using spline interpolation
time_smooth = np.linspace(time.min(), time.max(), 300)  # Create more points for smooth curve
spline = make_interp_spline(time, cell_concentration, k=3)  # k=3 for cubic smoothing
cell_concentration_smooth = spline(time_smooth)

# Identify Growth Phases
lag_phase_end = time[cell_concentration.diff().idxmax()] * 0.5  # Approximate Lag phase end
stationary_phase_start = time[cell_concentration.idxmax()]  # Start of stationary phase

# Plot Smoothed Cell Growth Curve
plt.figure(figsize=(8, 5))
sns.lineplot(x=time_smooth, y=cell_concentration_smooth, color='b', label="Cell Biomass Concentration", linewidth=2)
sns.scatterplot(x=time, y=cell_concentration, color='blue', label="Experimental Data")

# Highlight Phases
plt.axvspan(0, lag_phase_end, color='gray', alpha=0.3, label="Lag Phase")
plt.axvspan(lag_phase_end, stationary_phase_start, color='green', alpha=0.3, label="Exponential Phase")
plt.axvspan(stationary_phase_start, time.iloc[-5], color='blue', alpha=0.3, label="Stationary Phase")
plt.axvspan(time.iloc[-5], time.iloc[-1], color='red', alpha=0.3, label="Falling Phase")

# Labels & Aesthetics
plt.xlabel("Time (hrs)")
plt.ylabel("Cell Growth Concentration (mg/mL)")
plt.title("Growth Phases of Cell Biomass Over Time")
plt.legend()
plt.savefig(f"{folder}/cell_growth_curve (phases).png", dpi=300, bbox_inches="tight")
# plt.show()

# Print Estimated Phase Boundaries
print(f"Lag Phase ends around: {lag_phase_end:.2f} hrs")
print(f"Stationary Phase starts around: {stationary_phase_start:.2f} hrs")





# Bioflocculant Production - Experimental Data
bioflocculant_data = {
    "Time (hrs)": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0],
    "Absorbance": [0.48, 0.61, 0.75, 0.87, 0.96, 1.00, 1.03, 1.05, 1.07, 0.95, 0.84, 0.77, 0.66]
}

# Convert to DataFrame
bioflocculant_df = pd.DataFrame(bioflocculant_data)

# Convert Absorbance to Bioflocculant Concentration using Calibration Curve
bioflocculant_df["Bioflocculant Concentration (mg/mL)"] = (bioflocculant_df["Absorbance"] - intercept) / slope

# Compute Production Rate (dp/dt) using Finite Differences
bioflocculant_df["Production Rate (mg/mL per hr)"] = np.gradient(
    bioflocculant_df["Bioflocculant Concentration (mg/mL)"], bioflocculant_df["Time (hrs)"]
)

# Generate Smooth Time Points for Interpolation
smooth_time = np.linspace(bioflocculant_df["Time (hrs)"].min(), bioflocculant_df["Time (hrs)"].max(), 200)

# Interpolate for smooth curves
spline_concentration = make_interp_spline(bioflocculant_df["Time (hrs)"], bioflocculant_df["Bioflocculant Concentration (mg/mL)"])
smooth_concentration = spline_concentration(smooth_time)

spline_production_rate = make_interp_spline(bioflocculant_df["Time (hrs)"], bioflocculant_df["Production Rate (mg/mL per hr)"])
smooth_production_rate = spline_production_rate(smooth_time)

# Plot Bioflocculant Concentration Over Time (Smooth)
plt.figure(figsize=(8, 5))
sns.lineplot(x=smooth_time, y=smooth_concentration, color='b', label="Bioflocculant Concentration", linewidth=2)
sns.scatterplot(x=bioflocculant_df["Time (hrs)"], y=bioflocculant_df["Bioflocculant Concentration (mg/mL)"], color='b', label="Observed Data")
plt.xlabel("Time (hrs)")
plt.ylabel("Bioflocculant Concentration (mg/mL)")
plt.title("Bioflocculant Concentration Over Time")
plt.legend()
plt.savefig(f"{folder}/bioflocculant_concentration_curve.png", dpi=300)
# plt.show()

# Plot Bioflocculant Production Rate Over Time (Smooth)
plt.figure(figsize=(8, 5))
sns.lineplot(x=smooth_time, y=smooth_production_rate, color='green', label="Production Rate (dp/dt)", linewidth=2)
sns.scatterplot(x=bioflocculant_df["Time (hrs)"], y=bioflocculant_df["Production Rate (mg/mL per hr)"], color='green', label="Observed Data")
plt.xlabel("Time (hrs)")
plt.ylabel("Production Rate (mg/mL per hr)")
plt.title("Bioflocculant Production Rate Over Time")
plt.legend()
plt.savefig(f"{folder}/bioflocculant_production_rate.png", dpi=300)
# plt.show()

# Save DataFrame to CSV
bioflocculant_df.to_csv("bioflocculant_production_data.csv", index=False)
print("Data and plots saved successfully!")

# Luedeking-Piret Model
# Extract absorbance related to cell growth (OD600)
absorbance = cell_growth_df["Absorbance"].iloc[1:].values  # Ignore the first time point

# Compute (μ × Biomass)
mu_X_values = mu_values[1:] * absorbance  # Element-wise multiplication

# Extract production rate from bioflocculant data
production_rate = bioflocculant_df["Production Rate (mg/mL per hr)"].iloc[1:].values

# Fit Linear Model (Luedeking-Piret Equation: rp = α (μX) + β)
slope, intercept, r_value, _, _ = linregress(mu_X_values, production_rate)

# Predictions for plotting (linear model)
mu_X_pred = np.linspace(min(mu_X_values), max(mu_X_values), 100)
production_rate_pred = slope * mu_X_pred + intercept  # Linear equation

plt.figure(figsize=(10, 6))
sns.scatterplot(x=mu_X_values, y=production_rate, label="Observed Data", color='r')
sns.lineplot(x=mu_X_pred, y=production_rate_pred, label="Fitted Model", linestyle="--", color="b")
plt.xlabel("Specific Growth Rate × Biomass")
plt.ylabel("Production Rate (mg/mL per hr)")
plt.title("Luedeking-Piret Model Fit")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.7)

# 📌 Equation Annotation (Non-Distracting)
equation_text = (
    f"$r_p = {slope:.4f} (\\mu X) + {intercept:.4f}$\n"
    f"$R^2 = {r_value**2:.4f}$"
)

plt.annotate(
    equation_text,
    xy=(0.05, 0.85),  # Adjust position (relative to axes)
    xycoords="axes fraction",
    fontsize=12, color="black",
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='black')
)

# Save the plot
plt.savefig(f"{folder}/luedeking_piret_model.png", dpi=300)
plt.tight_layout()
# plt.show()

# Display estimated parameters
print(f"Estimated α (Growth-Associated Term): {slope:.4f}")
print(f"Estimated β (Non-Growth-Associated Term): {intercept:.4f}")
print(f"R-squared value: {r_value**2:.4f}")

print(f"Luedeking-Piret Model Equation: rp = {slope:.4f} * (μX) + {intercept:.4f}")



# Calibration Curve Data (Glucose)
calibration_data = {
    "Absorbance": [0, 0.08, 0.15, 0.32, 0.61, 0.96],  # Absorbance Values
    "Concentration (mg/mL)": [0, 0.2, 0.4, 0.6, 0.8, 1]  # Concentration Values
}

# Convert to DataFrame
sub_calib_df = pd.DataFrame(calibration_data)

# Perform Linear Regression (y = mx + b)
slope, intercept, r_value, _, _ = linregress(sub_calib_df["Concentration (mg/mL)"], sub_calib_df["Absorbance"])

# Print Calibration Equation
print(f"Calibration Curve Equation: Absorbance = ({slope:.4f} * Concentration) + {intercept:.4f}")

# Coefficient of determination
print(f"R-squared: {r_value**2:.4f}")

# Define x and y values
x = np.array(sub_calib_df["Concentration (mg/mL)"])
y = np.array(sub_calib_df["Absorbance"])

# Generate Smooth Curve (Cubic Spline)
x_smooth = np.linspace(x.min(), x.max(), 200)  # More points for smooth curve
spline = make_interp_spline(x, y, k=3)  # k=3 for cubic spline
y_smooth = spline(x_smooth)

# --- Plot 1: Scatter + Linear Fit ---
plt.figure(figsize=(8, 5))
plt.scatter(x, y, color='blue', label="Calibration Data")

# Add R-squared value to the graph
r_squared = r_value ** 2
plt.text(0.1, 0.8, f"R² = {r_squared:.4f}", fontsize=12, color="black", transform=plt.gca().transAxes)

plt.plot(x, slope * x + intercept, color='red', linestyle='--', label=f"Fit: y={slope:.4f}x + {intercept:.4f}")
plt.xlabel("Concentration (mg/mL)")
plt.ylabel("Absorbance")
plt.title("Calibration Curve: Absorbance vs. Concentration (Glucose)")
plt.legend()
plt.grid()
plt.savefig(f"{folder}/calibration_linear_glucose.png", dpi=300, bbox_inches="tight")  # Save Image
# plt.show()

# --- Plot 2: Scatter + Smooth Curve ---
plt.figure(figsize=(8, 5))
plt.scatter(x, y, color='blue', label="Calibration Data")

# Add R-squared value to the graph
r_squared = r_value ** 2
plt.text(0.1, 0.8, f"R² = {r_squared:.4f}", fontsize=12, color="black", transform=plt.gca().transAxes)

plt.plot(x_smooth, y_smooth, color='black', linestyle='-', linewidth=2, label="Cubic Spline Fit")
plt.plot(x, slope * x + intercept, color='red', linestyle='--', label=f"Fit: y={slope:.4f}x + {intercept:.4f}")
plt.xlabel("Concentration (mg/mL)")
plt.ylabel("Absorbance")
plt.title("Calibration Curve: Absorbance vs. Concentration (Glucose)")
plt.legend()
plt.grid()
plt.savefig(f"{folder}/calibration_smooth_glucose.png", dpi=300, bbox_inches="tight")  # Save Image
# plt.show()



# Substrate consumption 
# Experimental Data (time, absorbance, pH, temp)
experimental_data = {  
    "Time (hrs)": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0],  
    "Absorbance": [0.62, 0.61, 0.60, 0.59, 0.57, 0.53, 0.47, 0.41, 0.36, 0.51, 0.61, 0.73, 0.86]  
}  

substrate_df = pd.DataFrame(experimental_data)

# Convert Absorbance to Substrate Consumption Concentration using Calibration curve 
substrate_df["Substrate Concentration (mg/mL)"] = (substrate_df["Absorbance"] - intercept) / slope

substrate_df['Substrate Consumption Rate (mg/mL per hr)'] = -np.gradient(substrate_df["Substrate Concentration (mg/mL)"], substrate_df["Time (hrs)"])

substrate_df.head(2)

# Extract time and substrate concentration
time = substrate_df["Time (hrs)"].values
substrate_concentration = substrate_df["Substrate Concentration (mg/mL)"].values
consumption_rate = substrate_df["Substrate Consumption Rate (mg/mL per hr)"].values

# Create smooth curves using spline interpolation
time_smooth = np.linspace(time.min(), time.max(), 200)  # Generate finer time points
substrate_spline = make_interp_spline(time, substrate_concentration, k=3)  # Cubic spline
consumption_spline = make_interp_spline(time, consumption_rate, k=3)

substrate_smooth = substrate_spline(time_smooth)
consumption_smooth = consumption_spline(time_smooth)

# Plot Substrate Concentration Over Time
plt.figure(figsize=(8, 5))
sns.lineplot(x=time_smooth, y=substrate_smooth, label="Substrate Concentration", color="b", linewidth=2)
sns.scatterplot(x=time, y=substrate_concentration, color='b', edgecolor='black', label="Data Point")
plt.xlabel("Time (hrs)")
plt.ylabel("Substrate Concentration (mg/mL)")
plt.title("Substrate Concentration Over Time")
plt.legend()
plt.savefig(f"{folder}/substrate_concentration_over_time.png", dpi=300)
# plt.show()

# Plot Substrate Consumption Rate Over Time
plt.figure(figsize=(8, 5))
sns.lineplot(x=time_smooth, y=consumption_smooth, label="Consumption Rate", color="green", linewidth=2)
sns.scatterplot(x=time, y=consumption_rate, color='green', edgecolor='black', label="Data Point")
plt.xlabel("Time (hrs)")
plt.ylabel("Substrate Consumption Rate (mg/mL per hr)")
plt.title("Substrate Consumption Rate Over Time")
plt.legend()
plt.savefig(f"{folder}/substrate_consumption_rate.png", dpi=300)
# plt.show()

# Define the Monod model equation
def monod_model(s, umax, Ks):
    """
    Monod Equation: μ = (μmax * s) / (Ks + s)
    
    Parameters:
    S  = Substrate Concentration (mg/mL)
    umax = Maximum Specific Growth Rate
    Ks = Half-Saturation Constant
    
    Returns:
    Growth rate as per Monod model.
    """
    return (umax * s) / (Ks + s)  

# Fit the Monod model
popt, _ = curve_fit(monod_model, substrate_df["Time (hrs)"], substrate_df["Substrate Concentration (mg/mL)"])
umax, Ks = popt  # Extract fitted parameters

# Generate smooth curve for the model
t_smooth = np.linspace(min(substrate_df["Time (hrs)"]), max(substrate_df["Time (hrs)"]), 100)
predicted_concentration = monod_model(t_smooth, umax, Ks)

# Plot Experimental Data and Monod Model Fit
plt.figure(figsize=(8, 5))
sns.scatterplot(x=substrate_df["Time (hrs)"], y=substrate_df["Substrate Concentration (mg/mL)"], color='blue', label="Experimental Data")
sns.lineplot(x=t_smooth, y=predicted_concentration, color='red', linestyle='--', label="Monod Model Fit")

# Add equation and parameters as text
equation_text = f"Monod Model:\nμ = (umax * S) / (Ks + S)\n\nEstimated Parameters:\n- umax = {umax:.4f}\n- Ks = {Ks:.4f}"
plt.text(0.6 * max(substrate_df["Time (hrs)"]), 0.8 * max(substrate_df["Substrate Concentration (mg/mL)"]), 
         equation_text, fontsize=10, bbox=dict(facecolor='white', alpha=0.6))

# Labels and title
plt.xlabel("Time (hrs)")
plt.ylabel("Substrate Concentration (mg/mL)")
plt.title("Monod Model Fitting (Substrate Concentration vs Time)")
plt.legend()
plt.savefig(f"{folder}/monod_model_fit.png", dpi=300)
# plt.show()

substrate_df.to_csv('substrate_data.csv', index=False)