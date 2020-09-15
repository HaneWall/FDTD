import numpy as np

try:
    import torch
    torch.set_default_dtype (torch.float64)
    torch_available = True
    cuda_available = torch.cuda.is_available()
except ImportError:
    torch_available = False
    cuda_available = False

class Backend:
    #which bd is active?
    def __repr__(self):
        return self.__class__.__name__

class NumpyBackend(Backend):
    int = np.int64
    float = np.float64

    sin = staticmethod(np.sin)
    cos = staticmethod(np.cos)
    exp = staticmethod(np.exp)
    cosh = staticmethod(np.cosh)
    sum = staticmethod(np.sum)
    stack = staticmethod(np.stack)
    transpose = staticmethod(np.transpose)

    array = staticmethod(np.array)
    empty = staticmethod(np.empty)
    ones = staticmethod(np.ones)
    zeros = staticmethod(np.zeros)
    arange = staticmethod(np.arange)
    matmul = staticmethod(np.matmul)
    outer = staticmethod(np.outer)

if torch_available:
    class TorchBackend(Backend):
        int = torch.int64
        float = torch.float64

        sin = staticmethod(torch.sin)
        cos = staticmethod(torch.cos)
        exp = staticmethod(torch.exp)
        cosh = staticmethod(torch.cosh)
        sum = staticmethod(torch.sum)
        stack = staticmethod(torch.stack)
        transpose = staticmethod(torch.transpose)

        empty = staticmethod(torch.empty)
        ones = staticmethod(torch.ones)
        zeros = staticmethod(torch.zeros)
        arange = staticmethod(torch.arange)
        matmul = staticmethod(torch.matmul)
        outer = staticmethod(torch.ger)

        def array(self, arr, dtype = None):
            if dtype is None:
                dtype = torch.get_default_dtype()
                return torch.tensor(arr, device='cpu', dtype=dtype)
            return torch.is_tensor(arr)

        def convert_to_numpy(self, arr):
            if torch.is_tensor(arr):
                return arr.numpy()
            else:
                return np.asarray(arr)


#default:
backend = NumpyBackend()

def set_backend(name: str):
    '''if name == 'torch':
        torch.set_default_dtype(torch.float64)
        backend.__class__ = TorchBackend'''
    if name == 'numpy':
        backend.__class__ = NumpyBackend

    else:
        raise RuntimeError('backend?')