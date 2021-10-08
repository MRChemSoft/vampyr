from vampyr import LegendreBasis, InterpolatingBasis
import numpy as np
import matplotlib.pyplot as plt

class BasisPlotter:
    def __init__(self, basis_type="Interpolating", order=2):
        self.scaling_order = order

        if basis_type == "Interpolating":
            self.basis = InterpolatingBasis(self.scaling_order)
        if basis_type == "Legendre":
            self.basis = LegendreBasis(self.scaling_order)

    def plot_scaling(self, i=0, l=0, n=0):
        if i > self.scaling_order:
            print("Scaling order exceeded")
            return None
        x = np.arange(-0.5, 1.5, 0.001)
        phi = self.basis.scaling(i=i, l=l, n=n)

        phi_plt = [phi([x]) for x in x]
        return plt.plot(x, phi_plt)

    def plot_wavelet(self, i=0, l=0, n=0):
        if i > self.scaling_order:
            print("Scaling order exceeded")
            return None
        x = np.arange(-0.5, 1.5, 0.001)
        psi = self.basis.wavelet(i=i, l=l, n=n)

        psi_plt = [psi([x]) for x in x]
        return plt.plot(x, psi_plt)

if __name__ == "__main__":
    #Example of scaling plots

    basis = BasisPlotter(order=4)
    for i in range(5):
        basis.plot_scaling(i=i)
