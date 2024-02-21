import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the exponential function

def func(x, a, b):
    return a * x + b


# Generate sample data
x_data =np.array([2,3,4,5,6,7,8,9,10])
y_data = np.log([0.001088738459834888, 8.875312703468858e-05, 8.93275764168013e-06, 1.7630891112030712e-06, 
                 8.757034322377915e-06, 8.895232465026442e-07, 8.093380717715945e-08, 8.710266253689968e-08, 6.911504614261364e-09])
# Add some noise to the data
# np.random.seed(0)  # for reproducibility
# y_noise = 0.2 * np.random.normal(size=x_data.size)
# y_data += y_noise

# Fit the data to the exponential function
popt, pcov = curve_fit(func, x_data, y_data)

print(popt)

# Plot the original data and the fitted curve
plt.scatter(x_data, y_data, label='Data')
plt.plot(x_data, func(x_data, *popt), 'r-', label='Fitted curve')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Data Fitting')
plt.legend()
plt.grid(True)
plt.show()

# Print the parameters of the fitted exponential function
print("Fitted Parameters (a, b):", popt)
