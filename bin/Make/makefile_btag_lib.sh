if [[  $isSLC5 == "True" ]];
    then
    bash ${LQANALYZER_DIR}/bin/Make/make_btag_lib_c98.sh True
    
else
    bash ${LQANALYZER_DIR}/bin/Make/make_btag_lib_c11.sh True
fi