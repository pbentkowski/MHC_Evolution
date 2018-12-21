

# Stopka
cat << EOF
-------------------------------------------------------------------------------

End of calculations [$(date)].

-------------------------------------------------------------------------------
EOF

# Konczymy obliczenia, zawartosc katalogu $TMPDIR/output kopiujemy 
# do katalogu z ktorego zakolejkowano zadanie.
#mkdir $SLURM_SUBMIT_DIR/${OUTPUT_DIR}
cp -r $TMPDIR/* $SLURM_SUBMIT_DIR/

# Czyscimy katalog roboczy
rm -rf $TMPDIR
