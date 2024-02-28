import argparse
import os
import joblib
import numpy as np

def get_model_weights(model):
    weights = {}
    params = model.get_params()
    for param_name, param_value in params.items():
        if isinstance(param_value, (np.ndarray, list)):
            weights[param_name] = param_value
    return weights

def main(args):
    model_path = args.model
    output_dir = args.directory
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Load model
    model = joblib.load(model_path)
    print(model)

    # Get weights
    weights = get_model_weights(model)

    # Output weights
    for param_name, param_value in weights.items():
        print(f"Parameter: {param_name}")
        print(f"  Shape: {np.shape(param_value)}")
        # Saving weights to files
        weight_filename = f"{param_name}.npy"
        weight_path = os.path.join(output_dir, weight_filename)
        with open(weight_path, 'wb') as f:
            np.save(f, param_value)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract weights from a scikit-learn model.")
    parser.add_argument("model", type=str, help="Path to the scikit-learn model file.")
    parser.add_argument("directory", type=str, help="Directory to save weights.")
    args = parser.parse_args()
    main(args)
