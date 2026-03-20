"""
Implementation of molecular atlas mapping to standard brain template
Reference abagen.get_expression_data mapping strategy
"""

import os
import glob
import numpy as np
import nibabel as nib
from pathlib import Path
from typing import Union, Dict, List, Optional
import pandas as pd


class MolecularBrainMapper:
    """
    Molecular mapping tool class 
    
    Function:
        1. Load the molecular atlas data in .nii.gz format.
        2. Align to the standard brain template.
        3. Extract the statistical indicators at the ROI level.
        4. Handle missing values. 
    """
    
    def __init__(self, 
                 atlas_path: str,
                 data_dir: Optional[str] = None,
                 standard_template: str = 'MNI152',
                 verbose: bool = True):
        """
        Initialize mapper
        
        Parameters:
        ------
        atlas_path : str
            Standard brain template (atlas) file path, .nii format
        data_dir : str, optional
            Mitochondrial atlas data folder path
        standard_template : str
            Standard brain template type, such as 'MNI152', 'MNI305', etc.
        verbose : bool
            Whether to output detailed processing information
        """
        self.atlas_path = atlas_path
        self.data_dir = data_dir
        self.standard_template = standard_template
        self.verbose = verbose
        
        # Validate brain atlas file
        if not os.path.exists(atlas_path):
            raise FileNotFoundError(f"Brain atlas file does not exist: {atlas_path}")
        
        # Load brain atlas
        self.atlas = nib.load(atlas_path)
        self.atlas_data = self.atlas.get_fdata()
        self.roi_labels = np.unique(self.atlas_data[self.atlas_data > 0]).astype(int)
        
        if self.verbose:
            print(f"✓ Brain atlas loaded successfully")
            print(f"  - Atlas dimensions: {self.atlas_data.shape}")
            print(f"  - Number of ROIs: {len(self.roi_labels)}")
    
    def load_mitochondrial_data(self, 
                               file_pattern: str = "*.nii.gz") -> Dict[str, np.ndarray]:
        """
        Load molecular atlas files
        
        Parameters:
        ------
        file_pattern : str
            File matching pattern, default *.nii.gz
        
        Returns:
        ------
        dict
            {filename: nifti data array}
        """
        if self.data_dir is None:
            raise ValueError("data_dir not specified, please provide it during initialization")
        
        if not os.path.exists(self.data_dir):
            raise FileNotFoundError(f"Data directory does not exist: {self.data_dir}")
        
        mito_files = glob.glob(os.path.join(self.data_dir, file_pattern))
        
        if not mito_files:
            raise FileNotFoundError(f"No {file_pattern} files found in {self.data_dir}")
        
        mito_data = {}
        for file_path in sorted(mito_files):
            file_name = Path(file_path).stem.replace('.nii', '')
            try:
                img = nib.load(file_path)
                data = img.get_fdata()
                mito_data[file_name] = data
                
                if self.verbose:
                    print(f"✓ Loaded: {file_name} - Dimensions: {data.shape}")
            except Exception as e:
                print(f"✗ Failed to load: {file_name} - {str(e)}")
        
        return mito_data
    
    def map_to_atlas(self, 
                    mito_data: Dict[str, np.ndarray],
                    missing: str = 'centroids',
                    statistic: str = 'mean') -> pd.DataFrame:
        """
        Map molecular data to brain atlas (refer to abagen strategy)
        
        Parameters:
        ------
        mito_data : dict
            Molecular atlas data, {name: data array}
        missing : str
            Method to handle missing values
            - 'centroids': Use ROI centroids (default)
            - 'zeros': Fill with 0
            - 'nan': Keep as NaN
        statistic : str
            Extract ROI statistics
            - 'mean': Mean value
            - 'median': Median value
            - 'std': Standard deviation
            - 'max': Maximum value
            - 'min': Minimum value
        
        Returns:
        ------
        pd.DataFrame
            Mapping results at ROI level, rows are ROIs, columns are molecular types
        """
        if self.verbose:
            print(f"\nStart mapping data to brain atlas...")
            print(f"  - Missing value handling: {missing}")
            print(f"  - Statistic: {statistic}")
        
        # Initialize result matrix
        n_rois = len(self.roi_labels)
        result_matrix = np.zeros((n_rois, len(mito_data)))
        
        for col_idx, (file_name, data) in enumerate(mito_data.items()):
            # Check dimension consistency
            if data.shape != self.atlas_data.shape:
                print(f"Warning: {file_name} dimensions {data.shape} do not match brain atlas {self.atlas_data.shape}")
                data = self._resample_to_atlas(data)
            
            # Extract statistics for each ROI
            for roi_idx, roi_label in enumerate(self.roi_labels):
                roi_mask = self.atlas_data == roi_label
                roi_values = data[roi_mask]
                
                # Extract statistics
                if len(roi_values) > 0 and not np.all(np.isnan(roi_values)):
                    result_matrix[roi_idx, col_idx] = self._compute_statistic(
                        roi_values, statistic
                    )
                else:
                    # Handle missing values
                    result_matrix[roi_idx, col_idx] = self._handle_missing(
                        roi_label, roi_mask, data, missing
                    )
        
        # Convert to DataFrame
        df_result = pd.DataFrame(
            result_matrix,
            index=[f"ROI_{int(label)}" for label in self.roi_labels],
            columns=list(mito_data.keys())
        )
        
        if self.verbose:
            print(f"✓ Mapping completed")
            print(f"  - Result dimensions: {df_result.shape}")
            print(f"  - Number of ROIs: {len(df_result)}")
            print(f"  - Number of molecular types: {len(df_result.columns)}")
        
        return df_result
    
    def _resample_to_atlas(self, data: np.ndarray) -> np.ndarray:
        """
        Resample data to atlas space
        
        Parameters:
        ------
        data : np.ndarray
            Original data
        
        Returns:
        ------
        np.ndarray
            Resampled data
        """
        from scipy import ndimage
        
        # Use nearest-neighbor interpolation for scaling
        zoom_factors = np.array(self.atlas_data.shape) / np.array(data.shape)
        resampled = ndimage.zoom(data, zoom_factors, order=0)
        
        # Truncate to the same shape
        resampled = resampled[:self.atlas_data.shape[0],
                             :self.atlas_data.shape[1],
                             :self.atlas_data.shape[2]]
        
        return resampled
    
    @staticmethod
    def _handle_missing(roi_label: int, 
                       roi_mask: np.ndarray,
                       data: np.ndarray,
                       strategy: str) -> float:
        """
        Handle missing values
        
        Refer to abagen strategy:
        - 'centroids': Use the non-zero value closest to the ROI centroid
        - 'zeros': Fill with 0
        - 'nan': Keep as NaN
        """
        if strategy == 'zeros':
            return 0.0
        elif strategy == 'nan':
            return np.nan
        elif strategy == 'centroids':
            # Calculate ROI centroid
            coords = np.argwhere(roi_mask)
            if len(coords) > 0:
                centroid = coords.mean(axis=0).astype(int)
                centroid = np.clip(centroid, 0, np.array(data.shape) - 1)
                return data[tuple(centroid)]
            return np.nan
        else:
            return np.nan
        
    @staticmethod
    def _compute_statistic(values: np.ndarray, statistic: str) -> float:
        """Compute statistic"""
        values = values[~np.isnan(values)]  # Remove NaN values
        
        if len(values) == 0:
            return np.nan
        
        stats_dict = {
            'mean': np.mean,
            'median': np.median,
            'std': np.std,
            'max': np.max,
            'min': np.min,
            'sum': np.sum,
        }
        
        return stats_dict.get(statistic, np.mean)(values)
    def save_results(self, df_result: pd.DataFrame, output_path: str) -> None:
        """
        Save mapping results
        
        Parameters:
        ------
        df_result : pd.DataFrame
            Mapping result dataframe
        output_path : str
            Output file path
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if str(output_path).endswith('.csv'):
            df_result.to_csv(output_path)
        elif str(output_path).endswith('.xlsx'):
            df_result.to_excel(output_path)
        else:
            df_result.to_csv(output_path)
        
        if self.verbose:
            print(f"✓ Results saved to: {output_path}")