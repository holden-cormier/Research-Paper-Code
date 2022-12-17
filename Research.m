% Final Research Project Code - Holden Cormier




%% ONE MONTH
  % FIRST OPTION LOW TIME HORIZON AND STRIKE PRICE
clear all;
% ENTER MODEL PARAMETERS
  sigma = 0.77;                                          % annual volatility
  T     = 0.0833333;                                          % years to expiration
  K     = 180;                                          % option strike price
  p0    = 192.51;                                          % current asset price 
  r     = 0.0731;                                          % annual interest rate
  
% ENTER/COMPUTE TIME DISCRETIZATION PARAMETERS
  nt     = 30;
  deltat = T/(nt-1); 
  beta   = exp(-r*deltat);

% COMPUTE SHOCK DISTRIBUTION
 m = 5;                                                % number of shocks
 [e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];
  nx = length(x);

% PACK MODEL STRUCTURE
  clear model
  model.horizon = nt;                                   % time horizon
  model.func = 'pricing';                                % model functions DO I NEED TO CHANGE THIS?????
  model.discount = beta;                                % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions   
  model.discretestates = 2;                             % index of discrete states 
  model.params={K};                                     % other parameters

% DEFINE APPROXIMATION SPACE
  d    = 6;                                             % state halfrange in std deviations
  ns   = 500;                                           % degree of approximation
  s0   = log(p0);                                       % current log of asset price
  smin = s0 - d*sigma*sqrt(T);                          % compute minimum state
  smax = s0 + d*sigma*sqrt(T);                          % compute maximum state
  fspace = fundefn('lin',ns,smin,smax,[],[0;1]);        % function space  
  scoord = funnode(fspace);                             % state collocaton nodes
  s = gridmake(scoord); 
  
% INITIALIZE VALUE FUNCTION
  p = exp(s(:,1));                                      % compute prices at states
  v = zeros(2*ns,1);                                    % terminal value function

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  [~,s,v,x] = dpsolve(model,fspace,s,v);


% PLOT VALUE FUNCTION
  figure(1)
  p = exp(scoord{1});                               
  v = reshape(v,ns,2,nt+1); 
  v = squeeze(v(:,1,:));
  subplot(3,3,1)
  plot(p,v(:,1),p,max(p - K , 0)); 
  title('MongoDB Call Option Value Function (Strike $180, 1 Month)');
  xlabel('Asset Price'); ylabel('Value');



  
  
  
  
  
  
  
  
  
  
  
% SECOND OPTION LOW TIME HORIZON AND STRIKE PRICE
clear all;
% ENTER MODEL PARAMETERS
  sigma = 0.77;                                          % annual volatility
  T     = 0.0833333;                                          % years to expiration
  K     = 200;                                          % option strike price
  p0    = 192.51;                                          % current asset price 
  r     = 0.0731;                                          % annual interest rate
  
% ENTER/COMPUTE TIME DISCRETIZATION PARAMETERS
  nt     = 30;
  deltat = T/(nt-1);
  beta   = exp(-r*deltat);

% COMPUTE SHOCK DISTRIBUTION
 m = 5;                                                % number of shocks
 [e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];
  nx = length(x);

% PACK MODEL STRUCTURE
  clear model
  model.horizon = nt;                                   % time horizon
  model.func = 'pricing';                                % model functions DO I NEED TO CHANGE THIS?????
  model.discount = beta;                                % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions   
  model.discretestates = 2;                             % index of discrete states 
  model.params={K};                                     % other parameters

% DEFINE APPROXIMATION SPACE
  d    = 6;                                             % state halfrange in std deviations
  ns   = 500;                                           % degree of approximation
  s0   = log(p0);                                       % current log of asset price
  smin = s0 - d*sigma*sqrt(T);                          % compute minimum state
  smax = s0 + d*sigma*sqrt(T);                          % compute maximum state
  fspace = fundefn('lin',ns,smin,smax,[],[0;1]);        % function space  
  scoord = funnode(fspace);                             % state collocaton nodes
  s = gridmake(scoord); 
  
% INITIALIZE VALUE FUNCTION
  p = exp(s(:,1));                                      % compute prices at states
  v = zeros(2*ns,1);                                    % terminal value function

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  [c,s,v,x] = dpsolve(model,fspace,s,v);


% PLOT VALUE FUNCTION
  p = exp(scoord{1});                               
  v = reshape(v,ns,2,nt+1); 
  v = squeeze(v(:,1,:));
  subplot(3,3,2)
  plot(p,v(:,1),p,max(p - K , 0)); 
  title('MongoDB Call Option Value Function (Strike $200, 1 Month)');
  xlabel('Asset Price'); ylabel('Value');










% SECOND OPTION LOW TIME HORIZON AND STRIKE PRICE
clear all;
% ENTER MODEL PARAMETERS
  sigma = 0.77;                                          % annual volatility
  T     = 0.0833333;                                          % years to expiration
  K     = 260;                                          % option strike price
  p0    = 192.51;                                          % current asset price 
  r     = 0.0731;                                          % annual interest rate
  
% ENTER/COMPUTE TIME DISCRETIZATION PARAMETERS
  nt     = 30;
  deltat = T/(nt-1);
  beta   = exp(-r*deltat);

% COMPUTE SHOCK DISTRIBUTION
 m = 5;                                                % number of shocks
 [e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];
  nx = length(x);

% PACK MODEL STRUCTURE
  clear model
  model.horizon = nt;                                   % time horizon
  model.func = 'pricing';                                % model functions DO I NEED TO CHANGE THIS?????
  model.discount = beta;                                % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions   
  model.discretestates = 2;                             % index of discrete states 
  model.params={K};                                     % other parameters

% DEFINE APPROXIMATION SPACE
  d    = 6;                                             % state halfrange in std deviations
  ns   = 500;                                           % degree of approximation
  s0   = log(p0);                                       % current log of asset price
  smin = s0 - d*sigma*sqrt(T);                          % compute minimum state
  smax = s0 + d*sigma*sqrt(T);                          % compute maximum state
  fspace = fundefn('lin',ns,smin,smax,[],[0;1]);        % function space  
  scoord = funnode(fspace);                             % state collocaton nodes
  s = gridmake(scoord); 
  
% INITIALIZE VALUE FUNCTION
  p = exp(s(:,1));                                      % compute prices at states
  v = zeros(2*ns,1);                                    % terminal value function

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  [c,s,v,x] = dpsolve(model,fspace,s,v);


% PLOT VALUE FUNCTION
  p = exp(scoord{1});                               
  v = reshape(v,ns,2,nt+1); 
  v = squeeze(v(:,1,:));
  subplot(3,3,3)
  plot(p,v(:,1),p,max(p - K , 0)); 
  title('MongoDB Call Option Value Function (Strike $260, 1 Month) ');
  xlabel('Asset Price'); ylabel('Value');
















%%  SIX MONTHS
% FIRST OPTION SIX MONTH HORIZON AND $140 STRIKE
clear all;
% ENTER MODEL PARAMETERS
  sigma = 0.77;                                          % annual volatility
  T     = 0.5;                                          % years to expiration
  K     = 180;                                          % option strike price
  p0    = 192.51;                                          % current asset price 
  r     = 0.0731;                                          % annual interest rate
  
% ENTER/COMPUTE TIME DISCRETIZATION PARAMETERS
  nt     = 30;
  deltat = T/(nt-1);
  beta   = exp(-r*deltat);

% COMPUTE SHOCK DISTRIBUTION
 m = 5;                                                % number of shocks
 [e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];
  nx = length(x);

% PACK MODEL STRUCTURE
  clear model
  model.horizon = nt;                                   % time horizon
  model.func = 'pricing';                                % model functions DO I NEED TO CHANGE THIS?????
  model.discount = beta;                                % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions   
  model.discretestates = 2;                             % index of discrete states 
  model.params={K};                                     % other parameters

% DEFINE APPROXIMATION SPACE
  d    = 6;                                             % state halfrange in std deviations
  ns   = 500;                                           % degree of approximation
  s0   = log(p0);                                       % current log of asset price
  smin = s0 - d*sigma*sqrt(T);                          % compute minimum state
  smax = s0 + d*sigma*sqrt(T);                          % compute maximum state
  fspace = fundefn('lin',ns,smin,smax,[],[0;1]);        % function space  
  scoord = funnode(fspace);                             % state collocaton nodes
  s = gridmake(scoord); 
  
% INITIALIZE VALUE FUNCTION
  p = exp(s(:,1));                                      % compute prices at states
  v = zeros(2*ns,1);                                    % terminal value function

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  [c,s,v,x] = dpsolve(model,fspace,s,v);


% PLOT VALUE FUNCTION
  p = exp(scoord{1});                               
  v = reshape(v,ns,2,nt+1); 
  v = squeeze(v(:,1,:));
  subplot(3,3,4)
  plot(p,v(:,1),p,max(p - K , 0)); 
  title('MongoDB Call Option Value Function (Strike $180, 6 Months)');
  xlabel('Asset Price'); ylabel('Value');













  % SECOND OPTION HALF YEAR HORIZON AND STRIKE PRICE HORIZON AND STRIKE
  % PRICE &160
clear all;
% ENTER MODEL PARAMETERS
  sigma = 0.77;                                          % annual volatility
  T     = 0.5;                                          % years to expiration
  K     = 200;                                          % option strike price
  p0    = 192.51;                                          % current asset price 
  r     = 0.0731;                                          % annual interest rate
  
% ENTER/COMPUTE TIME DISCRETIZATION PARAMETERS
  nt     = 30;
  deltat = T/(nt-1);
  beta   = exp(-r*deltat);

% COMPUTE SHOCK DISTRIBUTION
 m = 5;                                                % number of shocks
 [e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];
  nx = length(x);

% PACK MODEL STRUCTURE
  clear model
  model.horizon = nt;                                   % time horizon
  model.func = 'pricing';                                % model functions DO I NEED TO CHANGE THIS?????
  model.discount = beta;                                % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions   
  model.discretestates = 2;                             % index of discrete states 
  model.params={K};                                     % other parameters

% DEFINE APPROXIMATION SPACE
  d    = 6;                                             % state halfrange in std deviations
  ns   = 500;                                           % degree of approximation
  s0   = log(p0);                                       % current log of asset price
  smin = s0 - d*sigma*sqrt(T);                          % compute minimum state
  smax = s0 + d*sigma*sqrt(T);                          % compute maximum state
  fspace = fundefn('lin',ns,smin,smax,[],[0;1]);        % function space  
  scoord = funnode(fspace);                             % state collocaton nodes
  s = gridmake(scoord); 
  
% INITIALIZE VALUE FUNCTION
  p = exp(s(:,1));                                      % compute prices at states
  v = zeros(2*ns,1);                                    % terminal value function

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  [c,s,v,x] = dpsolve(model,fspace,s,v);


% PLOT VALUE FUNCTION
  p = exp(scoord{1});                               
  v = reshape(v,ns,2,nt+1); 
  v = squeeze(v(:,1,:));
  subplot(3,3,5)
  plot(p,v(:,1),p,max(p - K , 0)); 
  title('MongoDB Call Option Value Function (Strike $200, 6 Months)');
  xlabel('Asset Price'); ylabel('Value');


























  % THIRD OPTION HALF YEAR HORIZON AND STRIKE PRICE $215
clear all;
% ENTER MODEL PARAMETERS
  sigma = 0.77;                                          % annual volatility
  T     = 0.5;                                          % years to expiration
  K     = 260;                                          % option strike price
  p0    = 192.51;                                          % current asset price 
  r     = 0.0731;                                          % annual interest rate
  
% ENTER/COMPUTE TIME DISCRETIZATION PARAMETERS
  nt     = 30;
  deltat = T/(nt-1);
  beta   = exp(-r*deltat);

% COMPUTE SHOCK DISTRIBUTION
 m = 5;                                                % number of shocks
 [e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];
  nx = length(x);

% PACK MODEL STRUCTURE
  clear model
  model.horizon = nt;                                   % time horizon
  model.func = 'pricing';                                % model functions DO I NEED TO CHANGE THIS?????
  model.discount = beta;                                % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions   
  model.discretestates = 2;                             % index of discrete states 
  model.params={K};                                     % other parameters

% DEFINE APPROXIMATION SPACE
  d    = 6;                                             % state halfrange in std deviations
  ns   = 500;                                           % degree of approximation
  s0   = log(p0);                                       % current log of asset price
  smin = s0 - d*sigma*sqrt(T);                          % compute minimum state
  smax = s0 + d*sigma*sqrt(T);                          % compute maximum state
  fspace = fundefn('lin',ns,smin,smax,[],[0;1]);        % function space  
  scoord = funnode(fspace);                             % state collocaton nodes
  s = gridmake(scoord); 
  
% INITIALIZE VALUE FUNCTION
  p = exp(s(:,1));                                      % compute prices at states
  v = zeros(2*ns,1);                                    % terminal value function

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  [c,s,v,x] = dpsolve(model,fspace,s,v);


% PLOT VALUE FUNCTION
  p = exp(scoord{1});                               
  v = reshape(v,ns,2,nt+1); 
  v = squeeze(v(:,1,:));
  subplot(3,3,6)
  plot(p,v(:,1),p,max(p - K , 0)); 
  title('MongoDB Call Option Value Function (Strike $260 , 6 Months) ');
  xlabel('Asset Price'); ylabel('Value');







































  %%  ONE YEAR
% FIRST OPTION ONE YEAR HORIZON AND $140 STRIKE
clear all;
% ENTER MODEL PARAMETERS
  sigma = 0.77;                                          % annual volatility
  T     = 1;                                          % years to expiration
  K     = 180;                                          % option strike price
  p0    = 192.51;                                          % current asset price 
  r     = 0.0731;                                          % annual interest rate
  
% ENTER/COMPUTE TIME DISCRETIZATION PARAMETERS
  nt     = 30;
  deltat = T/(nt-1);
  beta   = exp(-r*deltat);

% COMPUTE SHOCK DISTRIBUTION
 m = 5;                                                % number of shocks
 [e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];
  nx = length(x);

% PACK MODEL STRUCTURE
  clear model
  model.horizon = nt;                                   % time horizon
  model.func = 'pricing';                                % model functions DO I NEED TO CHANGE THIS?????
  model.discount = beta;                                % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions   
  model.discretestates = 2;                             % index of discrete states 
  model.params={K};                                     % other parameters

% DEFINE APPROXIMATION SPACE
  d    = 6;                                             % state halfrange in std deviations
  ns   = 500;                                           % degree of approximation
  s0   = log(p0);                                       % current log of asset price
  smin = s0 - d*sigma*sqrt(T);                          % compute minimum state
  smax = s0 + d*sigma*sqrt(T);                          % compute maximum state
  fspace = fundefn('lin',ns,smin,smax,[],[0;1]);        % function space  
  scoord = funnode(fspace);                             % state collocaton nodes
  s = gridmake(scoord); 
  
% INITIALIZE VALUE FUNCTION
  p = exp(s(:,1));                                      % compute prices at states
  v = zeros(2*ns,1);                                    % terminal value function

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  [c,s,v,x] = dpsolve(model,fspace,s,v);


% PLOT VALUE FUNCTION
  p = exp(scoord{1});                               
  v = reshape(v,ns,2,nt+1); 
  v = squeeze(v(:,1,:));
  subplot(3,3,7)
  plot(p,v(:,1),p,max(p - K , 0)); 
  title('MongoDB Call Option Value Function (Strike $180, 1 Year)');
  xlabel('Asset Price'); ylabel('Value');













  % SECOND OPTION ONE YEAR HORIZON AND STRIKE PRICE HORIZON AND STRIKE
  % PRICE &160
clear all;
% ENTER MODEL PARAMETERS
  sigma = 0.77;                                          % annual volatility
  T     = 1;                                          % years to expiration
  K     = 200;                                          % option strike price
  p0    = 192.51;                                          % current asset price 
  r     = 0.0731;                                          % annual interest rate
  
% ENTER/COMPUTE TIME DISCRETIZATION PARAMETERS
  nt     = 30;
  deltat = T/(nt-1);
  beta   = exp(-r*deltat);

% COMPUTE SHOCK DISTRIBUTION
 m = 5;                                                % number of shocks
 [e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];
  nx = length(x);

% PACK MODEL STRUCTURE
  clear model
  model.horizon = nt;                                   % time horizon
  model.func = 'pricing';                                % model functions DO I NEED TO CHANGE THIS?????
  model.discount = beta;                                % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions   
  model.discretestates = 2;                             % index of discrete states 
  model.params={K};                                     % other parameters

% DEFINE APPROXIMATION SPACE
  d    = 6;                                             % state halfrange in std deviations
  ns   = 500;                                           % degree of approximation
  s0   = log(p0);                                       % current log of asset price
  smin = s0 - d*sigma*sqrt(T);                          % compute minimum state
  smax = s0 + d*sigma*sqrt(T);                          % compute maximum state
  fspace = fundefn('lin',ns,smin,smax,[],[0;1]);        % function space  
  scoord = funnode(fspace);                             % state collocaton nodes
  s = gridmake(scoord); 
  
% INITIALIZE VALUE FUNCTION
  p = exp(s(:,1));                                      % compute prices at states
  v = zeros(2*ns,1);                                    % terminal value function

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  [c,s,v,x] = dpsolve(model,fspace,s,v);


% PLOT VALUE FUNCTION
  p = exp(scoord{1});                               
  v = reshape(v,ns,2,nt+1); 
  v = squeeze(v(:,1,:));
  subplot(3,3,8)
  plot(p,v(:,1),p,max(p - K , 0)); 
  title('MongoDB Call Option Value Function (Strike $200, 1 Year)');
  xlabel('Asset Price'); ylabel('Value');


























  % THIRD OPTION ONE YEAR HORIZON AND STRIKE PRICE $215
clear all;
% ENTER MODEL PARAMETERS
  sigma = 0.77;                                          % annual volatility
  T     = 1;                                          % years to expiration
  K     = 260;                                          % option strike price
  p0    = 192.51;                                          % current asset price 
  r     = 0.0731;                                          % annual interest rate
  
% ENTER/COMPUTE TIME DISCRETIZATION PARAMETERS
  nt     = 30;
  deltat = T/(nt-1);
  beta   = exp(-r*deltat);

% COMPUTE SHOCK DISTRIBUTION
 m = 5;                                                % number of shocks
 [e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];
  nx = length(x);

% PACK MODEL STRUCTURE
  clear model
  model.horizon = nt;                                   % time horizon
  model.func = 'pricing';                                % model functions DO I NEED TO CHANGE THIS?????
  model.discount = beta;                                % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions   
  model.discretestates = 2;                             % index of discrete states 
  model.params={K};                                     % other parameters

% DEFINE APPROXIMATION SPACE
  d    = 6;                                             % state halfrange in std deviations
  ns   = 500;                                           % degree of approximation
  s0   = log(p0);                                       % current log of asset price
  smin = s0 - d*sigma*sqrt(T);                          % compute minimum state
  smax = s0 + d*sigma*sqrt(T);                          % compute maximum state
  fspace = fundefn('lin',ns,smin,smax,[],[0;1]);        % function space  
  scoord = funnode(fspace);                             % state collocaton nodes
  s = gridmake(scoord); 
  
% INITIALIZE VALUE FUNCTION
  p = exp(s(:,1));                                      % compute prices at states
  v = zeros(2*ns,1);                                    % terminal value function

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  [c,s,v,x] = dpsolve(model,fspace,s,v);


% PLOT VALUE FUNCTION
  p = exp(scoord{1});                               
  v = reshape(v,ns,2,nt+1); 
  v = squeeze(v(:,1,:));
  subplot(3,3,9)
  plot(p,v(:,1),p,max(p - K , 0)); 
  title('MongoDB Call Option Value Function (S trike $260 , 1 Year) ');
  xlabel('Asset Price'); ylabel('Value');



