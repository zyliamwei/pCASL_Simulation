clear all
clc
close all

%% Part 1: Simulation of the effect of artery radius and flow velocity on labeling efficiency
% Inversion efficiency (IE) table from the original MRM paper reporting
% BIC-pCASL MRI; Magn Reson Med 2025; DOI: 10.1002/mrm.70086.
inveff=[0.5	0.3855
0.6	0.4220
0.7	0.4549
0.8	0.4848
0.9	0.5121
1.0	0.5370
1.1	0.5597
1.2	0.5807
1.3	0.5999
1.4	0.6174
1.5	0.6334
1.6	0.6483
1.7	0.6622
1.8	0.6751
1.9	0.6871
2.0	0.6983
2.1	0.7085
2.2	0.7183
2.3	0.7274
2.4	0.7359
2.5	0.7441
2.6	0.7515
2.7	0.7586
2.8	0.7653
2.9	0.7717
3.0	0.7777
3.1	0.7837
3.2	0.7888
3.3	0.7940
3.4	0.7992
3.5	0.8040
3.6	0.8083
3.7	0.8126
3.8	0.8170
3.9	0.8209
4.0	0.8248
4.1	0.8281
4.2	0.8315
4.3	0.8349
4.4	0.8384
4.5	0.8413
4.6	0.8442
4.7	0.8471
4.8	0.8500
4.9	0.8524
5.0	0.8553
5.1	0.8577
5.2	0.8601
5.3	0.8625
5.4	0.8643
5.5	0.8666
5.6	0.8684
5.7	0.8709
5.8	0.8726
5.9	0.8744
6.0	0.8761
6.1	0.8778
6.2	0.8795
6.3	0.8806
6.4	0.8823
6.5	0.8839
6.6	0.8852
6.7	0.8861
6.8	0.8877
6.9	0.8888
7.0	0.8901
7.1	0.8908
7.2	0.8920
7.3	0.8932
7.4	0.8938
7.5	0.8951
7.6	0.8961
7.7	0.8967
7.8	0.8981
7.9	0.8982
8.0	0.8997
8.1	0.9005
8.2	0.9012
8.3	0.9020
8.4	0.9027
8.5	0.9035
8.6	0.9043
8.7	0.9050
8.8	0.9060
8.9	0.9060
9.0	0.9077
9.1	0.9077
9.2	0.9083
9.3	0.9096
9.4	0.9101
9.5	0.9113
9.6	0.9115
9.7	0.9120
9.8	0.9133
9.9	0.9136
10.0	0.9154
10.1	0.9159
10.2	0.9163
10.3	0.9177
10.4	0.9183
10.5	0.9187
10.6	0.9202
10.7	0.9202
10.8	0.9209
10.9	0.9213
11.0	0.9228
11.1	0.9237
11.2	0.9242
11.3	0.9256
11.4	0.9267
11.5	0.9267
11.6	0.9272
11.7	0.9287
11.8	0.9298
11.9	0.9298
12.0	0.9304
12.1	0.9319
12.2	0.9331
12.3	0.9331
12.4	0.9338
12.5	0.9352
12.6	0.9352
12.7	0.9365
12.8	0.9373
12.9	0.9373
13.0	0.9386
13.1	0.9399
13.2	0.9399
13.3	0.9407
13.4	0.9407
13.5	0.9420
13.6	0.9432
13.7	0.9432
13.8	0.9441
13.9	0.9441
14.0	0.9453
14.1	0.9453
14.2	0.9463
14.3	0.9463
14.4	0.9472
14.5	0.9472
14.6	0.9482
14.7	0.9482
14.8	0.9490
14.9	0.9490
15.0	0.9499
15.1	0.9499
15.2	0.9507
15.3	0.9507
15.4	0.9511
15.5	0.9511
15.6	0.9518
15.7	0.9518
15.8	0.9518
15.9	0.9525
16.0	0.9525
16.1	0.9524
16.2	0.9524
16.3	0.9528
16.4	0.9528
16.5	0.9528
16.6	0.9532
16.7	0.9532
16.8	0.9532
16.9	0.9524
17.0	0.9524
17.1	0.9524
17.2	0.9524
17.3	0.9524
17.4	0.9524
17.5	0.9524
17.6	0.9524
17.7	0.9508
17.8	0.9508
17.9	0.9508
18.0	0.9502
18.1	0.9502
18.2	0.9502
18.3	0.9497
18.4	0.9497
18.5	0.9497
18.6	0.9471
18.7	0.9471
18.8	0.9471
18.9	0.9471
19.0	0.9457
19.1	0.9457
19.2	0.9457
19.3	0.9446
19.4	0.9446
19.5	0.9446
19.6	0.9446
19.7	0.9408
19.8	0.9408
19.9	0.9408
20.0	0.9383
20.1	0.9383
20.2	0.9383
20.3	0.9383
20.4	0.9364
20.5	0.9364
20.6	0.9364
20.7	0.9364
20.8	0.9311
20.9	0.9311
21.0	0.9311
21.1	0.9311
21.2	0.9273
21.3	0.9273
21.4	0.9273
21.5	0.9273
21.6	0.9273
21.7	0.9244
21.8	0.9244
21.9	0.9244
22.0	0.9244
22.1	0.9174
22.2	0.9174
22.3	0.9174
22.4	0.9174
22.5	0.9174
22.6	0.9120
22.7	0.9120
22.8	0.9120
22.9	0.9120
23.0	0.9120
23.1	0.9079
23.2	0.9079
23.3	0.9079
23.4	0.9079
23.5	0.9079
23.6	0.8989
23.7	0.8989
23.8	0.8989
23.9	0.8989
24.0	0.8989
24.1	0.8915
24.2	0.8915
24.3	0.8915
24.4	0.8915
24.5	0.8915
24.6	0.8915
24.7	0.8860
24.8	0.8860
24.9	0.8860
25.0	0.8860
25.1	0.8860
25.2	0.8860
25.3	0.8746
25.4	0.8746
25.5	0.8746
25.6	0.8746
25.7	0.8746
25.8	0.8746
25.9	0.8650
26.0	0.8650
26.1	0.8650
26.2	0.8650
26.3	0.8650
26.4	0.8650
26.5	0.8578
26.6	0.8578
26.7	0.8578
26.8	0.8578
26.9	0.8578
27.0	0.8578
27.1	0.8578
27.2	0.8437
27.3	0.8437
27.4	0.8437
27.5	0.8437
27.6	0.8437
27.7	0.8437
27.8	0.8437
27.9	0.8315
28.0	0.8315
28.1	0.8315
28.2	0.8315
28.3	0.8315
28.4	0.8315
28.5	0.8315
28.6	0.8315
28.7	0.8222
28.8	0.8222
28.9	0.8222
29.0	0.8222
29.1	0.8222
29.2	0.8222
29.3	0.8222
29.4	0.8222
29.5	0.8049
29.6	0.8049
29.7	0.8049
29.8	0.8049
29.9	0.8049
30.0	0.8049
];

% Testing a series of arterial radius and velocity
mi=1;
for v0=3:0.1:30             % Flow velocity from 3 to 30 cm/s
    wIE=[];
    wvel=[];
    for R=100:2:600         % arterial radius; from 200 to 600 um at steps of 10 um
        r=0:1:R;            % location within the artery by reference to the arterial center;
                            % r = 0 represents the center and r = R represents the
                            % arterial wall
        vel=v0.*(1-r.^2./R.^2);     % Laminar flow
        pos=find(vel<0.5);      % Matching the labeling efficiency table
        r(pos)=[];
        vel(pos)=[];
        IE=inveff(fix((vel-0.4)./0.1),2).';
        area=2*pi*r*1;       % Simualtion spatial resolution is 1 um
        wIE=[wIE sum(IE.*area)./sum(area)];
        wvel=[wvel sum(vel.*area)./sum(area)]; 
    end
    gwIE(mi,:)=wIE;
    mi=mi+1;
end
R=100:10:600;
% Plot the labeling efficiency map
figure;
imshow(gwIE,[0.6 1],'initialmag','fit');
colormap(hot);colorbar;

% Examine the statistical effect of artery radius and flow velocity
s_IE=gwIE(:);
s_rad=ones(size(gwIE,1),1)*(100:2:600);
s_rad=s_rad(:)*1e-6;
s_vel=(3:0.1:30).'*ones(1,size(gwIE,2));
s_vel=s_vel(:);
ds=dataset(s_IE,s_rad,s_vel);
lm1=fitlm(ds,'s_IE ~ s_rad + s_vel')
coefCI(lm1)

%% Part 2: Simulation of the cross-sectional profile of velocity distribution and IE
R=200;      % Artery radius of 200 um as an example
vmax=15;    % Flow velocity of 15 cm/s as an example
velmap=zeros(500,500);
velmap2=zeros(500,500);
for mi=1:500
    for ni=1:500
        if (mi-250).^2+(ni-250).^2 < R.^2 
            % Laminar flow 
            velmap(mi,ni)=(1-((mi-250).^2+(ni-250).^2)./R.^2).*vmax;
            % Plug flow
            velmap2(mi,ni)=vmax;
        end
    end
end
% Display the velocity maps for laminar and plug flows
figure;imshow(velmap,[0 20],'initialmag','fit');colormap(jet);colorbar;title('Velocity map: Laminar');
figure;imshow(velmap2,[0 20],'initialmag','fit');colormap(jet);colorbar;title('Velocity map: Plug');
% Converting the velocity matrix for IE by reference to table
iemap=zeros(500,500);
iemap2=zeros(500,500);
for mi=1:500
    for ni=1:500
        if velmap(mi,ni)>3      % Cutoff velocity for simulation
            % Laminar flow
            iemap(mi,ni)=inveff(fix((velmap(mi,ni)-0.4)./0.1),2);
            % Plug flow
            iemap2(mi,ni)=inveff(fix((velmap2(mi,ni)-0.4)./0.1),2);
        end
    end
end
% Display the IE maps for laminar and plug flows
figure;imshow(iemap,[0.6 1],'initialmag','fit');colormap(hot);colorbar;title('Labeling efficiency: Laminar');
figure;imshow(iemap2,[0.6 1],'initialmag','fit');colormap(hot);colorbar;title('Labeling efficiency: Plug');


%% Part 3: Simulation of dependence of labeling efficiency on flow patterns
R=200;                  % Artery radius of 200 um as an example
wIE=[];
for v0=3:0.1:30         % Velocity range of 3 to 30 cm/s
    r=0:1:R;            % location within the artery by reference to the arterial center;
                        % r = 0 represents the center and r = R represents the
                        % arterial wall
    vel=v0.*(1-r.^2./R.^2);     % Laminar flow
    pos=find(vel<0.5);           % Matching the range of flow velocity in table
    r(pos)=[];
    vel(pos)=[];
    IE=inveff(fix((vel-0.4)./0.1),2).';
    area=2*pi*r*1;       % the simualtion resolution is 1 um
    wIE=[wIE sum(IE.*area)./sum(area)];
end
% Display the dependence of IE on flow velocity for laminar flow
v0=3:0.1:30;
figure;plot(v0,wIE,'b-');
xlabel('Peak velocity (cm/s)');
ylabel('Averaged IE');
W_Plot;
ylim([0.5 1]);
% Display the dependence of IE on flow velocity for plug flow
hold on;plot(v0,inveff(fix((v0-0.4)./0.1),2).','r-');
legend('Laminar flow','Plug flow');


