'''
Pyfftw executor that is used in all 3D FFTs.

At the moment, batched FFT transforms of 3 3D arrays are performed.
'''
import os

import pyfftw

from quvac import config


class FFTExecutor:
    def __init__(self, grid_shape, nthreads=None, fft_axes=(0,1,2)):
        self.grid_shape = grid_shape
        self.nthreads = nthreads if nthreads else os.cpu_count()
        self.fft_axes = fft_axes

    def allocate_fft(self):
        if getattr(self, "tmp_shape", None) is None:
            self._allocate_fft()

    def _allocate_fft(self):
        # self.tmp_shape = (3,) + self.grid_shape
        self.tmp_shape = self.grid_shape
        self.tmp = pyfftw.zeros_aligned(self.tmp_shape, dtype=config.CDTYPE)
        
        self.forward_fftw = pyfftw.FFTW(
                self.tmp,
                self.tmp,
                axes=self.fft_axes,
                direction="FFTW_FORWARD",
                flags=(config.FFTW_FLAG,),
                threads=self.nthreads,
        )

        self.backward_fftw = pyfftw.FFTW(
                self.tmp,
                self.tmp,
                axes=self.fft_axes,
                direction="FFTW_BACKWARD",
                flags=(config.FFTW_FLAG,),
                threads=self.nthreads,
        )


def setup_fftw_executor(fft_executor, grid_shape):
    if fft_executor is None:
        fft_executor = FFTExecutor(grid_shape)
    fft_executor.allocate_fft()
    return fft_executor
