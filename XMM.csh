#! /bin/tcsh
#####################
# This scripts downloads all the XMM data that we use for this project.
# Some of the data is still propriotary, so replace _set_password_here
# below with the actual password.
#
# In step two, it performs the usual steps of data reprocessing with
# odfingest, cifbuild, and appropriate standard chains.
# For pointed observations, we use xmmextractor, which processes the
# proposal target. In some cases we don't need that so it's a waste of time
# but it's still convenient to treat all targets the same here and just call
# a single script with default parameters.
#
# The total run time for all the reduction is more than 2 days, so this script
# could be broken into parts to make some of the processing run in parallel.
#
# Call this script in the base data directory and initialize SAS before.
######################

# I found that the following system packages required by SAS were
# missing from my system. I just record them here so I don't have to
# look up again which perl error message corresponds to which CentOS
# package.
#
# sudo yum install perl-Switch
# sudo yum install perl-XML-LibXML

# XMM
foreach obsid ( 0670380501 0763880101 0201270101 141610601 0670380301 0670380601 0134540601 0134540701 0134540801 0134540901 0115830601 0115830701 0115840201 0115890601 0115890701 0115890901 0115900201 0115900601 0116150601 0116160601 0116160801 0116200601 0116200701 0116320801 0116340601 0116340701 0116710901 0116711001 0116890801 0116890901 0116891101 0791980201 0791980501 0129350201 0129350301 0134540101 0134540301 0134540401 0134540501 0117890801 0117890901 9307400002 9307400003 9077700007 9133200002 9159500004 9178700004 9214200006 9343900004 0044740601 0502360101 0502360201 0724690101 0801010301 0801010601 9116600004 9217000002 9300200003 0801010201 0784050101 0801010101 9117800002 9147300003 9208200004 9291600002 9310100002)
echo $obsid
curl -o $obsid.tar "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=$obsid&level=ODF"
mkdir $obsid
mkdir $obsid/odf
cd $obsid/odf
cp ../../$obsid.tar .
tar -xf $obsid.tar
tar -xf *.TAR
cd ../../
end

# Carl's program still needs username and password
foreach obsid ( 0843150101 0843150201 0843150301 0843150401 0843150501 0843151101 0843150601 0843150701 0843150801 0843150901 0843151001 )
  echo Downloading $obsid
  curl -o $obsid.tar "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=$obsid&AIOUSER=cmelis&AIOPWD=set_password_here&level=ODF"
  mkdir $obsid
  mkdir $obsid/odf
  cd $obsid/odf
  tar -xzf *.tar.gz
  tar -xf *.TAR
  cd ../../
end

# clean up tar files after unpacking saves > 20 GB of space!
find -name '*.TAR' -delete
find -name '*.tar.gz' -delete
rm *.tar

# process normal observations
# In first round, run for nominal source
foreach obsid ( `find . -maxdepth 1 -name '[0-8]*'` )
  echo processing $obsid
  cd $obsid
  setenv SAS_ODF ./odf
  setenv SAS_CCF ccf.cif
  xmmextractor >& screenoutput.log
  unsetenv SAS_ODF
  unsetenv SAS_CCF
  cd ..
end
  
# process slew
foreach obsid ( `find . -maxdepth 1 -name '9*'` )
  echo processing $obsid
  cd $obsid
  #setenv SAS_ODF ./odf
  #cifbuild
  setenv SAS_CCF ccf.cif
  #odfingest
  setenv SAS_ODF `ls -1 *SUM.SAS`
  #epproc
  setenv SAS_ATTITUDE RAF
  eslewchain
  unsetenv SAS_ATTITUDE
  unsetenv SAS_ODF
  unsetenv SAS_CCF
  cd ..
end
