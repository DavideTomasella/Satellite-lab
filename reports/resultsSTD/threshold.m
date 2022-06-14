clear all

peakFail = [0.0051, 0.0039, 0.0044, 0.0048];
averageFail = [4.0415e-08, 3.5564e-07, 8.4117e-08, 6.0324e-08];
meanSquareFail = [1.0275e-06, 1.0247e-06, 1.0256e-06, 1.0252e-06];
thSTDFail = (peakFail - averageFail)./sqrt(meanSquareFail-averageFail.^2)

peakSucc = [0.0223, 0.0133, 0.0109, 0.0098, 0.0066, 0.0051, 0.0724,...
            0.0453, 0.0283, 0.0187, 0.0117, 0.0084, 0.0086, 0.0055];
averageSucc = [2.4788e-07, 5.8117e-07, 4.4830e-07, 5.8958e-07, 2.0471e-07,...
                4.0415e-08, 9.5102e-08, 1.4174e-07, 2.2526e-07, 2.4671e-07,...
                2.1642e-07, 1.7954e-07, 1.7256e-07, 2.3481e-07];
meanSquareSucc = [1.0471e-06, 1.0348e-06, 1.0317e-06, 1.0318e-06, 1.0304e-06,...
                1.0275e-06, 1.3273e-06, 1.1449e-06, 1.0666e-06, 1.0425e-06,...
                1.0299e-06, 1.0238e-06, 1.0273e-06, 1.0267e-06]
thSTDSucc = (peakSucc - averageSucc)./sqrt(meanSquareSucc-averageSucc.^2)

figure;
histogram(round(thSTDSucc),'BinMethod','integers')
hold on;
histogram(round(thSTDFail),'BinMethod','integers');
grid on;
xticks((0:5:65));
legend("Successful tests","Failed tests");
title("Standard deviations threshold");
xlabel("Standard deviations");
ylabel("Number of occurrences");