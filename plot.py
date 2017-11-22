import matplotlib.pyplot as plt
from numpy import loadtxt

plt.style.use('ggplot')

f, axarr = plt.subplots(3, sharex=True, figsize=(6,10))

# Load all data
truth = loadtxt('truth.txt')
obs = loadtxt('obs.txt')
obs[obs == 0.0] = None
first_guess = loadtxt('first_guess.txt')
final_guess = loadtxt('final_guess.txt')

# Plot all data
for i in range(3):
    axarr[i].plot(truth[:,i], label="truth")
    axarr[i].plot(obs[:,i], 'x', label="observations", alpha=0.5)
    axarr[i].plot(first_guess[:,i], '--', label="first guess", alpha=0.5)
    axarr[i].plot(final_guess[:,i], ':', label="final guess")

axarr[2].legend()

plt.tight_layout()
plt.show()
