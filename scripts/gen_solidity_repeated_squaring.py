import json
import sys
import argparse


def generate_circom_code(n):
    # Check if n is a positive integer
    if n < 1:
        raise ValueError("The number of squarings 'n' must be a positive integer.")

    # Start with the template header
    circom_code = """/* Repeated squaring of x to get y */

template RepeatedSquaring() {

    signal input x;
    signal output y;

    /// repeated squaring of x to get y
"""

    # Initialize signals
    for i in range(n):
        circom_code += f"    signal tmp{i};\n"

    # Generate the squaring operations
    circom_code += f"    tmp0 <== x * x;\n"
    for i in range(1, n):
        circom_code += f"    tmp{i} <== tmp{i-1} * tmp{i-1};\n"

    # Set the output
    circom_code += f"    y <== tmp{n-1};\n"

    # End the template and main component
    circom_code += "}\n\ncomponent main = RepeatedSquaring();\n"

    # Save the circom code to a file
    file_name = f"circom/repeated_squaring_{n}.circom"
    with open(file_name, "w") as file:
        file.write(circom_code)

    print(f"The circom code has been saved to {file_name}")

    return circom_code


def generate_witness(n):
    # Calculate x and y values based on the number of squarings
    x = 3
    y = x
    for _ in range(n):
        y = y * y % 21888242871839275222246405745257275088548364400416034343698204186575808495617

    # Create the witness dictionary
    witness = {"x": str(x), "y": str(y)}

    # Save the witness to a JSON file
    file_name = f"circom/witness_{n}.json"
    with open(file_name, "w") as file:
        json.dump(witness, file, indent=4)

    print(f"The witness JSON has been saved to {file_name}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate circom code and witness for repeated squaring."
    )
    parser.add_argument("n", type=int, help="The number of squarings")
    args = parser.parse_args()

    n = args.n
    generate_circom_code(n)
    generate_witness(n)


if __name__ == "__main__":
    main()
