% MFDP04 Function file for American option pricing model
function out = pricing(flag,s,x,e,K);
switch flag
case 'f'; % REWARD FUNCTION
  out = (exp(s(:,1))-K).*x.*(1-s(:,2));
case 'g'; % STATE TRANSITION FUNCTION
  out = [s(:,1)+e  s(:,2)|x];
end