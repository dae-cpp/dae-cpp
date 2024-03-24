# dae-cpp

A simple but powerful header-only C++ solver for systems of Differential and Algebraic Equations (DAE).

## What is dae-cpp

### How does it work

### The main features of the solver

## Installation

## Testing

If you already cloned the project without `--recurse-submodules`, you can initialize and update `googletest` submodule by running

```bash
git submodule update --init
```

Then build and run the tests:

```bash
mkdir build
cd build
cmake ..
make
ctest
```

## How to use

### Step 0. Include dae-cpp into the project

### Step 1. Define the initial state vector

### Step 2. Set up the RHS

### Step 3. Set up the Mass matrix

### Step 4. Set up the Jacobian matrix

### Step 5. Set the solver options

### Step 6. Solve the system

#### Optional: Set up Observer

#### Optional: Set up Event function

### A note about Sparse Matrix Format

## Contribution and Feedback

Thank you for considering contributing to the project! Whether you're an experienced developer or just starting out, your ideas and improvements can make this project even better. No contribution is too small!

### How to contribute

0. Create a GitHub issue if you want to suggest or discuss your changes.
1. Fork the repository and clone it to your local machine.
2. Create a new branch for your contributions.
3. Make your changes and ensure they adhere to our coding standards.
4. Add tests and examples (if relevant), test your changes thoroughly.
5. Submit a pull request with a clear description of your changes and why they are beneficial.

### Feedback

Feel free to create a [GitHub issue](https://github.com/ikorotkin/dae-cpp/issues) for any questions, suggestions, or feedback you may have.

## Licensing

- [dae-cpp](https://github.com/ikorotkin/dae-cpp) is licensed under the [MIT license](https://github.com/ikorotkin/dae-cpp/blob/master/LICENSE).
- [autodiff](https://github.com/autodiff/autodiff) is licensed under the [MIT license](https://github.com/autodiff/autodiff/blob/main/LICENSE).
- [Eigen]() ...
