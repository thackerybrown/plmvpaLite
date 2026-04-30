% This make.m is for MATLAB and OCTAVE under Windows, Mac, and Unix
% this version is modified 12/12/25 by TIB with LDFLAGS="$LDFLAGS
% -ld_classic", which is a temporary fix known to help with mex compiling
% on current Matlab and Tahoe setup - not known when this will stop working

function make()
try
	% This part is for OCTAVE
	if(exist('OCTAVE_VERSION', 'builtin'))
		mex libsvmread.c
		mex libsvmwrite.c
		mex -I.. train.c linear_model_matlab.c ../linear.cpp ../newton.cpp ../blas/daxpy.c ../blas/ddot.c ../blas/dnrm2.c ../blas/dscal.c
		mex -I.. predict.c linear_model_matlab.c ../linear.cpp ../newton.cpp ../blas/daxpy.c ../blas/ddot.c ../blas/dnrm2.c ../blas/dscal.c
	% This part is for MATLAB
	% Add -largeArrayDims on 64-bit machines of MATLAB
	else
		mex LDFLAGS="$LDFLAGS -ld_classic" -largeArrayDims libsvmread.c
		mex LDFLAGS="$LDFLAGS -ld_classic" -largeArrayDims libsvmwrite.c
		mex LDFLAGS="$LDFLAGS -ld_classic" -I.. -largeArrayDims train.c linear_model_matlab.c ../linear.cpp ../newton.cpp ../blas/daxpy.c ../blas/ddot.c ../blas/dnrm2.c ../blas/dscal.c
		mex LDFLAGS="$LDFLAGS -ld_classic" -I.. -largeArrayDims predict.c linear_model_matlab.c ../linear.cpp ../newton.cpp ../blas/daxpy.c ../blas/ddot.c ../blas/dnrm2.c ../blas/dscal.c
	end
catch err
	fprintf('Error: %s failed (line %d)\n', err.stack(1).file, err.stack(1).line);
	disp(err.message);
	fprintf('=> Please check README for detailed instructions.\n');
end
