
%Change-point condition data
dataCP = learningAndPerception;

%Oddball condition data
dataOB = oddBall_Inference;

%% 

%Analysis of data

figure

plot([dataCP.X', dataOB.B],[dataCP.expectationX, dataOB.expectationB], '.')
ylabel('Model max-posterior value')
xlabel('ground truth value of stimulus')

figure 

plot([dataCP.shannon', dataOB.shannon], [dataCP.relative_entropy_Mu', dataOB.relative_entropy_C'], 'or', 'markerSize', 3);
xlabel('suprise of outcome, internal representation (Shannon information content)');
ylabel('Relative entropy last trial - current trial');

figure
subplot(2, 1, 1)
plot([-360,360], [-360,360], '--k');
plot([-360,360], [0,0], '--k');
plot([dataCP.predError, dataOB.predError], [dataCP.muUpdate,dataOB.cUpdate], 'or', 'markerSize', 6);
ylabel('update');
xlabel('prediction error');

xes = [ones(size([dataCP.predError, dataOB.predError]')), [dataCP.predError, dataOB.predError]', [dataCP.shannon; dataOB.shannon'] - mean([dataCP.shannon; dataOB.shannon']), [dataCP.predError, dataOB.predError]' .* ([dataCP.shannon; dataOB.shannon'] - mean([dataCP.shannon; dataOB.shannon']))];
[B,BINT,R,RINT,STATS]  = regress([dataCP.muUpdate,dataOB.cUpdate]', xes);
modelMu = sum(B' .* xes, 2);

subplot(2, 1, 2)
plot([-360,360], [-360,360], '--k');
plot([-360,360], [0,0], '--k');
plot(modelMu, [dataCP.muUpdate,dataOB.cUpdate],  'o', 'markerSize', 6, 'markerFaceColor', 'b');
% plot(modelMu, [dataCP.muUpdate,dataOB.cUpdate],  'o', 'markerSize', 6, 'markerFaceColor', 'b');
ylabel('Actual update')
xlabel('Regression-predicted update (predError, Shannon)')

figure;

plot([abs(dataCP.predError), abs(dataOB.predError)], [dataCP.shannon', dataOB.shannon], 'o', 'markerSize', 6, 'markerFaceColor', 'b' );
xlabel('Prediction error');
ylabel('Shannon information content');


figure; 

plot([dataCP.predictionErrorOnX', dataOB.predictionErrorOnB'], [dataCP.perceptualErrorOnX', dataOB.perceptualErrorOnB'], '.');
plot([-100, 100], [-100, 100], '--k')
plot([-100, 100], [0, 0], '--k')

xlabel('prediction error')
ylabel('perceptual error')