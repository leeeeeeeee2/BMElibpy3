#/bin/sh
#Uses the already created interfaces to build the .so files:

echo "Building the .so files to use the Python interface to the MVN library"

# Example of signature building
# f2py -m mvProAG2 mvProAG2.f dcuhre.f symInv.f only: mvproag2 -h mvProAG2.p
# Then edit the signature file and add the input parameters ( intent(in) ) and output parameters ( intent(out) )

#builds the mvnAG1 interface
echo "Compiling mvnAG1"
f2py2.3 mvnAG1.pyf -c mvnPack.f

# builds the uvProNR interface
echo "Compiling uvProNR"
f2py2.3 uvProNR.pyf -c uvProNR.f qromb3.f trapzd3.f polint3.f

# builds the uvMomNR interface
echo "Compiling uvMomNR"
f2py uvMomNR.pyf -c uvMomNR.f qromb3.f trapzd3.f polint3.f

# builds the uvMomVecAG2 interface
echo "Compiling uvMomVecAG2"
f2py mvMomVecAG2.pyf -c mvMomVecAG2.f dcuhre.f symInv.f

# builds the mvProAG2 interface
echo "Compiling mvProAG2"
f2py mvProAG2.pyf -c mvProAG2.f dcuhre.f symInv.f

echo "Mvnlib interfaces created"