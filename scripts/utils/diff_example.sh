for i in $(find . \! -type d | grep -v svn | grep -wv doc); do diff $i /data/kdi_prod/project_result/716/01.00/Hi-C_pipeline_test_ev_all2/$i || echo $i; done
