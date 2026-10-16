% QUICKSTART_FULL_CHAIN_OPTIMIZED Run and verify the optimized full FEC path.
root = fileparts(mfilename('fullpath'));
cd(root);
[result,report] = run_optimized_full_chain();
testReport = test_full_chain_throughput_optimization();
assert(testReport.passed && report.fec_valid_reconstruction_exact);
