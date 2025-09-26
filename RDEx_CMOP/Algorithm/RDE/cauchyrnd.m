function r = cauchyrnd(mu, gamma)
    r = mu + gamma * tan(pi * (rand - 0.5));
end