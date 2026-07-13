# SimulatingPixelResponse
The digitisation of the MightyPix (MP) requires dedicated event model classes that are capable of simulating the full detector response chain, from the initial energy deposit to the final detector readout and clustering stage.

## Fast-Track Reconstruction: 
- Implemented a numerically stable track approximation using polynomial functions.  
- Optimised performance by solving the chi2￼ minimisation via Cholesky Decomposition, exploiting matrix symmetry for speed. 
