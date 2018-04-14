rng(100);
% Defining network model parameters
vt = 1;                     % Spiking threshold
tau_m = 10;                 % Membrane time constant [ms]
g_m = 1;                    % Neuron conductance
Nsig = 0.25;                 % Variance amplitude of current
Nmean = 0.75;                % Mean current to neurons
tau_I = 10;                 % Time constant to filter the synaptic inputs
N = 500;                   % Number of neurons in total
NE = 0.5*N;                 % Number of excitatory neurons
NI = 0.5*N;                 % Number of inhibitory neurons
dt = 1;                     % Simulation time bin [ms]
T = 300/dt;                 % Simulation length 
W = 2/sqrt(N);                  % Connectivity strength

% Initialization

v = rand(N,1)*vt;           % membrane potential
vv = zeros(N,1);            % variable that notes if v crosses the threshold
Iback = zeros(N,1);         % building up the external current
SP = 0;                     % recording spike times
Ichem = zeros(N,1);         % current coming from synaptic inputs
Iext = zeros(N,1);     % external current
post_syn = zeros(T,1);
raster = [];                % save spike times for plotting

% loop over the time
for t = 1:T
    Iback = Iback + dt/tau_I*(-Iback +randn(N,1));          % generate a colored noise for the current
    Iext = Iback/sqrt(1/(2*(tau_I/dt)))*Nsig+Nmean;         % rescaling the noise current to have the correct mean and variance

    Ichem(1:NE) = Ichem(1:NE) + dt/tau_I*(-Ichem(1:NE) + W*(sum(vv(1:NE))-vv(1:NE))-W*(sum(vv(NE+1:end)))); % current to excitatory neurons coming from the synaptic inputs
    Ichem(NE+1:end) = Ichem(NE+1:end) + dt/tau_I*(-Ichem(NE+1:end) -W*(sum(vv(NE+1:end))-vv(NE+1:end))+W*(sum(vv(1:NE)))); % current to inhibitory neurons coming from the synaptic inputs
    Itot = Iext+Ichem;
    post_syn(t) = Ichem(50);
    %%%%%%%%%%% To insert integrate-and-fire model here  %%%%%%%%%%%%%
   
    v = v+ dt./tau_m.*(-v+g_m*Itot);         % Euler method for IF neuron
  %  if t == 200
   %     v(150:201) = vt;
   % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vv =(v>=vt);                                        % spike if voltage crosses the threshold    
    
    v = (1-vv).*v;                                      % reset after spike
    SP = find(vv);                                      % find the spike times
    raster=[raster;t*ones(length(SP),1),SP];            % save spike times for plotting
    %raster_1= [raster;t*ones(length(SP),1),SP];
end

raster_1 = raster
post_syn_1 = post_syn
% Plot the raster output
figure; hold on;
plot(1:300, post_syn_0,'r');
plot(1:300, post_syn_1,'g');
plot(1:300, post_syn_2,'m');
plot(1:300, post_syn_3,'b');
h = figure; hold on;
%plot(raster_0(:,1)*dt, raster_0(:,2),'.r')
%plot(raster_1(:,1)*dt, raster_1(:,2),'.b')
xlim([100 300])
xlabel('time [ms]','fontsize',20)
ylabel('neuron index','fontsize',20)
set(gca,'fontsize',20);
set(gca,'YDir','normal')