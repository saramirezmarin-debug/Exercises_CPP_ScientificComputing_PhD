T = readtable("../output/zipf.csv");

figure;
loglog(T.rank, T.frequency, "o-");
grid on;

xlabel("Rank");
ylabel("Frequency");

title("Zipf Law");