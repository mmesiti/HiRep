

MKFILE=../Make/MkFlags
cp $MKFILE $MKFILE.bk

sed -e 's/#.*$//;s/^\s*//' $MKFILE.bk | grep -vE '^\s*$|^REPR|^GAUGE_GROUP|^NG|^CFLAGS' > $MKFILE

cat $MKFILE

TEST_ALGEBRA_1(){
  GROUP=$1
  shift;
  for REPR in REPR_FUNDAMENTAL REPR_SYMMETRIC REPR_ANTISYMMETRIC REPR_ADJOINT
  do
    for NG in $@
    do
    if [ "$NG" -ne 2 ] || [ "$REPR" != "REPR_ANTISYMMETRIC" ]; then
      touch $MKFILE
      cd Algebra
      echo make REPR=$REPR GAUGE_GROUP=$GROUP NG=$NG CFLAGS="-Wall -std=c99 -g -O0" > $GROUP\_$REPR\_$NG\_check_algebra_1.out
      make REPR=$REPR GAUGE_GROUP=$GROUP NG=$NG CFLAGS="-Wall -std=c99 -g -O0" && \
      ./check_algebra_1 >> $GROUP\_$REPR\_$NG\_check_algebra_1.out
      cd ..

    fi
    done
  done
  
}

TEST_ALGEBRA_1 GAUGE_SUN 2 3 4 
TEST_ALGEBRA_1 GAUGE_SPN 2 4 6 8

cp $MKFILE.bk $MKFILE
