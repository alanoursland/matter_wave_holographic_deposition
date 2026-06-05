"""Torch device selection helper."""

import torch


def get_device() -> torch.device:
    """Return the best available torch device (CUDA if available, else CPU)."""
    return torch.device('cuda' if torch.cuda.is_available() else 'cpu')
