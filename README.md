# torontonian-julia
Julia implementation of the functions Torontonian and threshold detection probabilities (without displacement).

There are two implementations for the computation of the Torontonian function: the first one uses directly the definition, the second one uses the recursive algorithm introduced in [arxiv:2109.04528](https://arxiv.org/abs/2109.04528).

The recursive algorithm provides a polynomial speedup of the computation (when the number of modes is increased), as demonstrated in the following pictures:

![time_comparison](https://user-images.githubusercontent.com/95931266/190510160-44e99f6a-1bf6-4fcd-b1b4-b910b5a5f276.png)

![time_comparison_2](https://user-images.githubusercontent.com/95931266/190509883-9369e83a-d358-466e-b2fd-59f99453ef69.png)

To run the tests write `julia tor_tests.jl` on a terminal. 

