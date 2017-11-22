import matplotlib.pyplot as plt
from numpy import loadtxt

plt.style.use('ggplot')

f, axarr = plt.subplots(3, sharex=True, figsize=(6,10))

# Load all data
truth = loadtxt('truth.txt')
obs = loadtxt('obs.txt')
first_guess = loadtxt('first_guess.txt')
final_guess = loadtxt('final_guess.txt')

# Plot all data
for i in range(3):
    axarr[i].plot(truth[:,0], truth[:,i+1], label="truth")
    axarr[i].plot(obs[:,0], obs[:,i+1], 'x', label="observations", alpha=0.5)
    axarr[i].plot(first_guess[:,0], first_guess[:,i+1], '--', label="first guess", alpha=0.5)
    axarr[i].plot(final_guess[:,0], final_guess[:,i+1], ':', label="final guess")

axarr[2].legend()

plt.tight_layout()
plt.show()
