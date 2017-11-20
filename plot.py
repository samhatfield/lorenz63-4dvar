import matplotlib.pyplot as plt
from numpy import loadtxt

plt.style.use('ggplot')

# Plot truth
truth = loadtxt('truth.txt')
plt.plot(truth[:,0], label="truth")

# Plot observations
obs = loadtxt('obs.txt')
obs[obs == 0.0] = None
plt.plot(obs[:,0], 'x', label="observations", alpha=0.5)

# Plot first guess
truth = loadtxt('first_guess.txt')
plt.plot(truth[:,0], '--', label="first guess", alpha=0.5)

# Plot final guess
truth = loadtxt('final_guess.txt')
plt.plot(truth[:,0], ':', label="final guess")

plt.legend()

plt.show()
