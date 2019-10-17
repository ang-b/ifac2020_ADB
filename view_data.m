%% output
figure(1)
for i=1:N
    subplot(3,2,i);
    plot(t, reshape(yT(1,i,:), [1 nTot]), t, xref(1,:));
    ylim([-1 1]);
end

%% decoupled error
figure(2)
for i=1:N
    subplot(3,2,i);
    plot(t, reshape(yT(1,i,:) - xd(1,i,:), [1 nTot]));
end

%% coupled error
figure(3)
for i=1:N
    subplot(3,2,i);
    plot(t, reshape(yT(1,i,:) - xc(1,i,:), [1 nTot]));
end

%% residuals
figure(4)
for i=1:N
    subplot(3,2,i);
    plot(t, dcres(i,:).');
end

%% attacker
figure(5)
plot(t, xA)

%% inputs
figure(6)
for i=1:N
    subplot(3,2,i);
    plot(t, reshape(u(:,i,:), [mI nTot]));
end