import os
import glob
import pandas as pd
from skimage import io
from skimage.measure import regionprops_table

def process_image(mask_path, tiff_path, output_base):
    try:
        # Load the mask and Tiff image
        mask = io.imread(mask_path)
        tiff_image = io.imread(tiff_path)
        print(f"Processing mask: {mask_path}, TIFF: {tiff_path}")

        # Ensure dimensions match
        if tiff_image.shape[:2] != mask.shape:
            print(f"Error: Dimension mismatch between {mask_path} and {tiff_path}")
            return
        
        # Channel mapping for your specific needs (channels 1, 2, and 3)
        channel_mapping = {1: 'G6', 2: 'CD45_1', 3: 'CD45_2'}

        # Process channels 1, 2, and 3
        for i, name in channel_mapping.items():
            channel_image = tiff_image[..., i]

            # Measure properties using regionprops_table
            properties = regionprops_table(
                mask, 
                intensity_image=channel_image,
                properties=('label', 'area', 'mean_intensity')
            )
            output_path = f"{output_base}_{name}.csv"
            pd.DataFrame(properties).to_csv(output_path, index=False)

        print(f"Processed and saved for {tiff_path}")
    except Exception as e:
        print(f"Error processing {mask_path} and {tiff_path}: {e}")

def main():
    try:
        print("Script started.")
        working_directory = input("Enter the working directory: ").strip()
        tiff_files = glob.glob(os.path.join(working_directory, "*.tif"))
        print(f"Found {len(tiff_files)} TIFF files")

        if len(tiff_files) == 0:
            print("No TIFF files found.")
            return

        print("TIFF files found:")
        for tiff_path in tiff_files:
            print(tiff_path)
            base_name = os.path.splitext(os.path.basename(tiff_path))[0]
            mask_path = os.path.join(working_directory, f"{base_name}_cp_masks.png")
            print(f"Checking for mask: {mask_path}")

            if os.path.exists(mask_path):
                output_base = os.path.join(working_directory, base_name)
                print(f"Processing {tiff_path} with corresponding mask: {mask_path}")
                process_image(mask_path, tiff_path, output_base)
            else:
                print(f"Mask not found for {tiff_path}. Expected: {mask_path}")
    except Exception as e:
        print(f"Error in main function: {e}")

if __name__ == "__main__":
    main()