clear all
blue = [0 0.4470 0.7410];
red = [0.6350 0.0780 0.1840];
purple = [.63,.13,.94];
att = red;
rep = blue;
green = [0.4660 0.6740 0.1880];
syms x y z real
syms alpha beta gamma real
syms d1 d2 d3 real

d = [d1, d2, d3];
vf = [x*(alpha-alpha*x-(alpha+beta+d1)*y+d2*z); %[output:group:1b51c0be] %[output:8cbedc38]
    y*(beta+d1*x-beta*y-(beta+gamma+d3)*z); %[output:8cbedc38]
    z*(gamma-(gamma+alpha+d2)*x+d3*y-gamma*z)] %system (1.4) %[output:group:1b51c0be] %[output:8cbedc38]
F = [x;y;z;x+y+z-1] %[output:98e1ddc4]
K = F;
for i = 1:size(F,1)
    K(i) = collect(simplify(gradient(F(i),[x,y,z]).'*vf/F(i)),[x,y,z]);
end
K %[output:5032ddc0]


sp = solve(vf); sp = [sp.x,sp.y,sp.z] %[output:35b79825]

%%
num = [1,1,1,0,0,0] %parameters %[output:12b4514a]

vfi = subs(vf,[alpha,beta,gamma,d],num) %[output:1a9f4afb]

spi = solve(vfi);  spi = [spi.x,spi.y,spi.z] %equilibra %[output:8cda668a]

J = jacobian(vfi,[x,y,z]) %Jacobian matrix %[output:0341fa2d]
temp = [];
for i = 1:size(spi,1)
    if min(spi(i,:)) >= 0
        temp = [temp;i];
    end
end
spi = spi(temp,:) %[output:1d4f8e38]
lam = [];
ji = jacobian(vfi,[x,y,z]);
attP = []; repP = [];
for i = 1:size(spi,1) %[output:group:16c97af9]
    disp('------------------------------------------') %[output:52d4e928] %[output:60de2f77] %[output:2166523e] %[output:0666ac36] %[output:49bd0171]
    spi(i,:),jii = subs(ji,[x,y,z],spi(i,:)); %[output:2b7b0eed] %[output:8f519491] %[output:89d2ef0b] %[output:131321f4] %[output:41a89f9f]
    lam = double(expand(eig(jii))) %[output:63163dba] %[output:17d781c8] %[output:015dfcca] %[output:81a7e044] %[output:65ad05ee]
    if lam < 0
        attP = [attP,i];
    elseif lam > 0
        repP = [repP, i];
    end
end %[output:group:16c97af9]

vfi %[output:678614a5]
spi = double(spi) %[output:20369ba9]
p = spi(end,:);  %equilibrium on the invariant 2-simplex
%%
%[text] ## 
warning('off');
axisrange = [0,1.5,0,1.5,0,1.3];
figure; hold on; %[output:29cc3ebe]

space = linspace(0,1,100);
plot3(zeros(1,100),space*axisrange(4),zeros(1,100),'k--','linewidth',1.2); %[output:29cc3ebe]
plot3(zeros(1,100),zeros(1,100),space*axisrange(6),'k--','linewidth',1.2); %[output:29cc3ebe]
plot3(space*axisrange(2),zeros(1,100),zeros(1,100),'k--','linewidth',1.2); %[output:29cc3ebe]
text('Interpreter','latex','String','$O$','Position',[0,-0.2,-0.2], 'FontSize',20); %[output:29cc3ebe]
text('Interpreter','latex','String','$x_1$','Position',[axisrange(2)+.2,-.2,0.05],'FontSize',15); %[output:29cc3ebe]
text('Interpreter','latex','String','$x_2$','Position',[-.2,axisrange(4)-.01,0.14],'FontSize',15); %[output:29cc3ebe]
text('Interpreter','latex','String','$x_3$','Position',[0,-0.3,axisrange(6)],'FontSize',15); %[output:29cc3ebe]

%plot invariant 2-simplex
fill3([1,0,0],[0,1,0],[0,0,1],'--','FaceColor',[0.93,0.91,0.91],... %[output:29cc3ebe]
    'FaceAlpha',.6,'EdgeLighting',"none")%,'EdgeColor','none') %[output:29cc3ebe]

xi = linspace(0,1,100);
plot3(xi,1-xi,zeros(1,100),'k','LineWidth',1.2) %[output:29cc3ebe]
plot3(zeros(1,100),xi,1-xi,'k','LineWidth',1.2) %[output:29cc3ebe]
plot3(xi,zeros(1,100),1-xi,'k','LineWidth',1.2) %[output:29cc3ebe]


%plot equilibria
plot3(p(1),p(2),p(3),'o','markerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k') %[output:29cc3ebe]
plot3(1,0,0,'o','markerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k') %[output:29cc3ebe]
plot3(0,1,0,'o','markerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k') %[output:29cc3ebe]
plot3(0,0,1,'o','markerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k') %[output:29cc3ebe]
plot3(0,0,0,'o','markerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k') %[output:29cc3ebe]


% plot phase trajectories
p0 = [0.5,0.1,0.2];
p0(2) = 1-p0(1)-p0(3);
[~,Y] = ode45(@df,[0,20],p0); %[output:48a428ce] %[output:03c5b764] %[output:0eb999ef] %[output:4ca30baa] %[output:92245277] %[output:0bb6ff28] %[output:4789eab6] %[output:3ee060de] %[output:2f3f835c] %[output:1bd5d37c] %[output:29d8bd80] %[output:41994cd9] %[output:14e12ba0] %[output:4ef2c694] %[output:229b0ec2] %[output:3533ea85] %[output:3024aeae] %[output:778dab18] %[output:44077167] %[output:6c8c13ea] %[output:98431882] %[output:444e1e57] %[output:6a0a3cb0] %[output:695dbe4a] %[output:7e75c12c] %[output:3736a73d] %[output:8da9db11] %[output:5ce2ced0] %[output:29ab2444] %[output:6537de23] %[output:12c28644] %[output:1e32ade5] %[output:431d25eb] %[output:5ec6c8a5] %[output:5a243eb0] %[output:97d69e80] %[output:4b1b874a] %[output:3b4fcfde] %[output:53f1ecfa] %[output:3baaa464] %[output:9a9bac92] %[output:3229033f] %[output:2b4363d4] %[output:8fac1365] %[output:0e55d2bd] %[output:9e68206e] %[output:213f36bc] %[output:8aa8bcfd] %[output:7438ed62] %[output:72dc17bf] %[output:69b81093] %[output:228d8c36] %[output:0032b4f3] %[output:0c17e45c] %[output:25cba3f3] %[output:3fd518d5] %[output:5cb7dbe1] %[output:7ab91dce] %[output:9c3e44a8] %[output:3b9d0293] %[output:2482b4f8] %[output:0b5e2b34] %[output:23ce012b] %[output:7d477538] %[output:41a15f57] %[output:5849ee3d] %[output:8c00e294] %[output:12980bd2] %[output:920bc24b] %[output:47e00a74] %[output:76b0d150] %[output:298ab721] %[output:4671bddc]
plot3(Y(:,1),Y(:,2),Y(:,3),'color',red,'linewidth',1.2) %[output:29cc3ebe]

p0 = [4/3,4/3,5/3];
[~,Y] = ode45(@df,[0,30],p0); %[output:05438899] %[output:9daff058] %[output:1819c86d] %[output:69cd2354] %[output:7627c088] %[output:977f7072] %[output:26f20cde] %[output:6d106ede] %[output:2a36a45f] %[output:0dd4615f] %[output:1a197072] %[output:99ce0af3] %[output:3ad2a0b0] %[output:2652b34f] %[output:1b32b151] %[output:8a256509] %[output:5e31f127] %[output:187fc282] %[output:328340af] %[output:8a6c360c] %[output:32f384ef] %[output:79f830c2] %[output:0092d02a] %[output:2c043606] %[output:8bf0e7c2] %[output:2cc86898] %[output:6025bd09] %[output:90a0d0c5] %[output:499d8f8b] %[output:04674943] %[output:79786ad2] %[output:593b8dab] %[output:4dfb785b] %[output:0ea445d0] %[output:06f2268d] %[output:20bae2b7] %[output:8bffd7b6] %[output:44779c46] %[output:175c4043] %[output:1744e439] %[output:09310e8b] %[output:172a0f6e] %[output:542fd33f] %[output:7f8271f7] %[output:61558be3] %[output:443becdd] %[output:2ffaffeb] %[output:31834072] %[output:8da1575f] %[output:9b7e4d16] %[output:08c9c99e] %[output:91d75aa8] %[output:6d70e2d1] %[output:0837221b] %[output:26b4ebba] %[output:9fdb5700] %[output:9bed7a84] %[output:545ba9ce] %[output:7d721b90] %[output:3dc0f81c] %[output:57cabdfe] %[output:2809412a] %[output:6d1654a7] %[output:937e4eee] %[output:894e4de0] %[output:06076702] %[output:193ace1a] %[output:752c27a0] %[output:0bf87d0d] %[output:2ce7cd77] %[output:349ecc9e] %[output:6d837fed] %[output:7e03d500] %[output:0d70d53c] %[output:0ff8600d] %[output:3cf68df6] %[output:60c0e9e1] %[output:90a686fa] %[output:38d6e69d] %[output:43df9d4b] %[output:60eeefc8] %[output:52c04a16] %[output:58d167ec] %[output:8fcd971c] %[output:5c4a75e2] %[output:80d020c8] %[output:86766d06] %[output:8ed249d5] %[output:16d5abcd] %[output:254a2ed6] %[output:162d7526] %[output:052c4857] %[output:76bd546d] %[output:6e855607] %[output:75a5e231] %[output:69f8667d] %[output:91dba932] %[output:0f1f0578] %[output:728ee629] %[output:02a63255] %[output:52a7c19a] %[output:6bdead69] %[output:4a086848] %[output:296f398a] %[output:58d08c81] %[output:528ed5be] %[output:791f2c75] %[output:3cd0f9b7] %[output:3a0afff4] %[output:5173badd] %[output:288b760c] %[output:4d3db7be] %[output:9b884830] %[output:2c9e3764] %[output:4c88a3de] %[output:2a586522] %[output:22fc056a] %[output:8420ffbd] %[output:8746a4aa] %[output:5dff2e80] %[output:57f57f70] %[output:4f48c788] %[output:806adeed] %[output:78557b0c] %[output:3a79179a] %[output:3dcc7bd7] %[output:8f2a7996] %[output:3a76b020] %[output:07e4bb0a] %[output:409d7227] %[output:43a981e4] %[output:967da8f4] %[output:71dd68fc]
plot3(Y(:,1),Y(:,2),Y(:,3),'color',red,'linewidth',1.2) %[output:29cc3ebe]

p0 = [4/3,3/4,1];
[~,Y] = ode45(@df,[0,30],p0); %[output:54dcedb5] %[output:3729288e] %[output:3cdfe3ff] %[output:4c60964d] %[output:73921f30] %[output:98eb2766] %[output:181d86fd] %[output:5db6f4b3] %[output:38bdcca9] %[output:5facba9e] %[output:4491e213] %[output:69e11642] %[output:115999ff] %[output:8a331cf7] %[output:657675d8] %[output:0bcaf2f0] %[output:2800d1f3] %[output:22b1eaef] %[output:4e9a7cfc] %[output:748a19e1] %[output:6fee66d3] %[output:1f7e7fd3] %[output:98aa73e6] %[output:50010a4b] %[output:11b4cf06] %[output:7eb4fb62] %[output:30724f9c] %[output:9a959871] %[output:5ec4f0aa] %[output:79e9f0ec] %[output:838b18cc] %[output:85273423] %[output:4867c0ce] %[output:6b1c2fe6] %[output:4071a987] %[output:2e73f8b1] %[output:59e9086a] %[output:5610f497] %[output:69fb9678] %[output:5ee42cad] %[output:696dccc7] %[output:6319a07a] %[output:2007325d] %[output:03bd9cec] %[output:1805b338] %[output:91b5a663] %[output:3ee16d6e] %[output:632b57bf] %[output:5c839f59] %[output:9b396e31] %[output:658f5540] %[output:1a73b7a9] %[output:157b6084] %[output:873ecacb] %[output:4807c517] %[output:27bdb534] %[output:8d746b0f] %[output:6c0db045] %[output:996777e0] %[output:5a4154b0] %[output:7415f0e7] %[output:862abcbf] %[output:057a0996] %[output:5432e3f9] %[output:11622cfd] %[output:6bd84909] %[output:174f43ed] %[output:40e848aa] %[output:7a6c7b5e] %[output:1be63df6] %[output:3bc0e1ce] %[output:0d01b568] %[output:72acb21e] %[output:77afd7a1] %[output:2e617915] %[output:9887abbc] %[output:73aa143c] %[output:3e1a8a95] %[output:6fd3cf39] %[output:65f764ff] %[output:3cff0c2b] %[output:3105e8c6] %[output:63ff485e] %[output:67247625] %[output:320dc3e4] %[output:0b3c89f4] %[output:4790f2c4] %[output:7163a3c1] %[output:297b940c] %[output:029a9f03] %[output:81578a73] %[output:9068849d] %[output:40e226e7] %[output:5b20dded] %[output:8769e9ad] %[output:80023480] %[output:15b51bd2] %[output:21c5a378] %[output:6510b820] %[output:168f3912] %[output:04e959d6] %[output:54fc490b] %[output:0a917969] %[output:17692bf1] %[output:57170966] %[output:83157faf] %[output:664da028] %[output:3b07c81c] %[output:3d09fde4] %[output:1eb1ef8e] %[output:613ec346] %[output:538ee531] %[output:259d9abe] %[output:096c9d0c] %[output:322d5ac4] %[output:62fe6b0a] %[output:8a774e59] %[output:3e4f288b] %[output:7750f67d] %[output:743cee12] %[output:61fdd38d] %[output:78ccc79d] %[output:031a10f9] %[output:075c14d5] %[output:331b52b2] %[output:271daa25] %[output:2345ac56] %[output:0ca1c8e9] %[output:492118ce] %[output:6e84c9b1] %[output:855260fc] %[output:5b585b0c] %[output:0229dd8c] %[output:6c041a82] %[output:9993e8bb] %[output:8f7afa11] %[output:86e82cc5] %[output:86cca90c] %[output:58d3010e]
plot3(Y(:,1),Y(:,2),Y(:,3),'color',red,'linewidth',1.2) %[output:29cc3ebe]

p0 = [0.1,0.1,0.11];
[~,Y] = ode45(@df,[0,30],p0); %[output:2faab95e] %[output:2e645e16] %[output:68c0d57e] %[output:70a18445] %[output:178883dd] %[output:9831720f] %[output:54e4a787] %[output:996ba7af] %[output:5e41a2e7] %[output:00ce9e37] %[output:0b8e9600] %[output:840fa490] %[output:9773f19e] %[output:8d613736] %[output:81572fc1] %[output:56784ed5] %[output:575b9f6f] %[output:9bbfc1e0] %[output:573c6ea0] %[output:1c7d4558] %[output:0930483b] %[output:1157a434] %[output:3c18aeab] %[output:5c5814ee] %[output:7964d2fe] %[output:97b9b513] %[output:111287c5] %[output:8d7d6ed9] %[output:860c331b] %[output:25f2f9cc] %[output:054c0ea1] %[output:2e87fcbb] %[output:489ab95c] %[output:4a595c2b] %[output:964f0080] %[output:2eae40fc] %[output:7292aaa6] %[output:2d9effa8] %[output:32c37d9f] %[output:27a486b0] %[output:77ccb615] %[output:71186c69] %[output:2b0b8a30] %[output:32086bb9] %[output:29277dbd] %[output:3502ea4a] %[output:74b7d7ac] %[output:61adffe2] %[output:5f3fdf94] %[output:7f42369c] %[output:05f9d35a] %[output:38f90cc0] %[output:721db9dd] %[output:3a965c1f] %[output:6325deb8] %[output:067775cd] %[output:103f627b] %[output:6c5d2ae3] %[output:85013286] %[output:9f509151] %[output:8bcbb9f7] %[output:3e1dc70d] %[output:396d5105] %[output:31c8de88] %[output:75d3f1c7] %[output:0d5aa878] %[output:1996b230] %[output:82e06869] %[output:62fdfd67] %[output:1c3a9219] %[output:9f5269c2] %[output:19803bc1] %[output:2563ab05] %[output:4ab734e4] %[output:86d8f36b] %[output:830bc724] %[output:94bc8491] %[output:99cbf1d5] %[output:03316099] %[output:26cb6bfe] %[output:9896a64d] %[output:886524af] %[output:3db67aab] %[output:85e54546] %[output:337c6343] %[output:0a549124] %[output:24969e3c] %[output:366fdbc5] %[output:96ba7690] %[output:32688f23] %[output:7a174383] %[output:3a3b4447] %[output:6f060927] %[output:12f89498] %[output:717d6dc1] %[output:7adede2f] %[output:98319ac6] %[output:03f75bb0] %[output:9a4527e4] %[output:7abe1f36] %[output:958d5cf2] %[output:9d16d211] %[output:9fed4249]
plot3(Y(:,1),Y(:,2),Y(:,3),'color',blue,'linewidth',1.2) %[output:29cc3ebe]
[~,Y] = ode45(@df_pos,[0,40],p0); %[output:393dab2d] %[output:69aafaf9] %[output:0dfb7868] %[output:56aaa37b] %[output:59d5ca29] %[output:3e97accf] %[output:2d4449c5] %[output:368e07aa] %[output:6feaeb46] %[output:85ba7983] %[output:49cd605d] %[output:0a4d4382] %[output:008394c8] %[output:4b8f7a17] %[output:630aa356] %[output:323fb556] %[output:6855e9ad] %[output:6c59756f] %[output:69570c30] %[output:7f2dd004] %[output:60edf16e] %[output:7b7980cd] %[output:9fa91884] %[output:5fe36a5c] %[output:2357fd62] %[output:8bdefd2b] %[output:5d5b39ff] %[output:5fdd333c] %[output:67cfe10b] %[output:961911c1] %[output:5fbcfd8a] %[output:4ee6508d] %[output:71cfed9d] %[output:6288a73a] %[output:9ba075ef] %[output:3be775b8] %[output:47d0543c] %[output:0352f9f2] %[output:98c62e4f] %[output:35c6f0e9] %[output:3e42848f] %[output:17bbf2c5] %[output:890b6cac] %[output:68bb5fa3] %[output:73a1a15d] %[output:4f3edcf2] %[output:31726c99] %[output:63bc83fd] %[output:6b62f62e] %[output:468bc253] %[output:996af791] %[output:01b4315b] %[output:87524361] %[output:7dc7fe0e] %[output:758368d4] %[output:2a336769] %[output:5ae0ff06] %[output:0f572356] %[output:7b038bf6] %[output:1a0ab645] %[output:603ae915] %[output:3fe04ca0] %[output:751bd134] %[output:8cedb0b8] %[output:1ca336f1] %[output:3a52a283] %[output:447875c5] %[output:70223ab6] %[output:349c4c29] %[output:58f2c191] %[output:48100f33] %[output:560398d4] %[output:676dcdce] %[output:62549ce2] %[output:5789374a] %[output:29a86abe] %[output:014f258f] %[output:537dba83] %[output:687e1cac] %[output:605e4f11] %[output:82fcb9a2] %[output:164a1b9b] %[output:8680788e] %[output:9d007e56] %[output:1411562e] %[output:4c9a98db] %[output:99053e30] %[output:94b32cc8] %[output:2e2989b8] %[output:6d7360bb] %[output:3e163704] %[output:6299b67e] %[output:30dce336] %[output:3977cbf6] %[output:91984ce2] %[output:85e9b295] %[output:8395638a] %[output:4b7c1319] %[output:38912115] %[output:8d275658] %[output:04dc5795] %[output:28f8d4eb] %[output:95761ab5] %[output:2f32a4be] %[output:0e01745a] %[output:52c7210e] %[output:0adfd303] %[output:8976b4ba] %[output:20e85172] %[output:4c42c122] %[output:1c2ceb15] %[output:5c8b0d75] %[output:6efe9c59] %[output:215dd0f7] %[output:1500daaa] %[output:774cbab4] %[output:1acfea74] %[output:3742a461] %[output:22d724dd] %[output:0c96a232] %[output:0f9564f2] %[output:64cc28c1] %[output:80d7af2f] %[output:6acc7d40] %[output:99178dba] %[output:9baf54c6] %[output:6cbe9aa7] %[output:4fec3391] %[output:370ea9d6] %[output:96c6f664] %[output:36a74d55] %[output:63540d10] %[output:4a3572f7] %[output:31ecf745] %[output:17a1ca9d] %[output:9dc2287a] %[output:026938f3] %[output:23bb3469] %[output:17b7f131]
plot3(Y(:,1),Y(:,2),Y(:,3),'color',blue,'linewidth',1.2) %[output:29cc3ebe]

p0 = [0.05,0.2,0.05];
[~,Y] = ode45(@df,[0,30],p0); %[output:66a36a1f] %[output:780f318d] %[output:2e891c64] %[output:73b549b2] %[output:2ccb067b] %[output:9a259f81] %[output:72f55d9b] %[output:89e24c96] %[output:3cb5ab61] %[output:87457d59] %[output:49fdd6af] %[output:603e8026] %[output:6750cbbc] %[output:9c2e7b92] %[output:155dba2e] %[output:91509444] %[output:19b84c56] %[output:95ced293] %[output:9163499f] %[output:2eb91684] %[output:57edf667] %[output:0d7079cf] %[output:10cda76a] %[output:3bd85b4a] %[output:37cd5761] %[output:4f56fa4d] %[output:49de5903] %[output:8ab9444d] %[output:51c5d869] %[output:14038319] %[output:3b5bd3f1] %[output:297df7ba] %[output:615b3336] %[output:1bb6b617] %[output:3522fa7b] %[output:27033efb] %[output:8fc6eb8b] %[output:91b0ffca] %[output:97f567d4] %[output:55e615f0] %[output:8c45b033] %[output:93ce9353] %[output:29a973bf] %[output:0eaf5d2d] %[output:56a241f9] %[output:40463efd] %[output:4962730c] %[output:9ca36a46] %[output:1930a20f] %[output:77dbc6ec] %[output:4e78437a] %[output:6255467b] %[output:9d7a9c01] %[output:239da291] %[output:835a22e7] %[output:3feb62ff] %[output:8d7c810b] %[output:146b3e42] %[output:7a8346e7] %[output:05a0b632] %[output:56f865bc] %[output:8f5b2dac] %[output:68416f78] %[output:151afdef] %[output:55fd7c1b] %[output:73177822] %[output:730f8e2e] %[output:94e58b25] %[output:5fe03690] %[output:49c96240] %[output:08fd31bd] %[output:9b33ae04] %[output:06bccb80] %[output:18521f56] %[output:00718c4b] %[output:8de07be5] %[output:949214ed] %[output:57b07787] %[output:24140659] %[output:33eaa987] %[output:97090215] %[output:68fd236b] %[output:228b0167] %[output:2b338af4] %[output:245d99ad] %[output:50193d57] %[output:530c57ba] %[output:747100ed] %[output:77d57786] %[output:21471de4] %[output:6a97357f] %[output:287a8c89] %[output:022a3a31] %[output:0e47ca2f] %[output:24546a6f] %[output:2e02dce4] %[output:2b1188e9] %[output:1b4dee9e] %[output:95f7baee] %[output:0a43a5d7] %[output:84c26120] %[output:5c19dc26] %[output:9c8b6f0a] %[output:3d50784f] %[output:1f31e725] %[output:433f0e23] %[output:2c975871] %[output:433a8922] %[output:30b3ae83] %[output:94f26261] %[output:5cc1fc40] %[output:28c31db4] %[output:17fc6145] %[output:6aba45e6] %[output:87e6c965] %[output:99e14172] %[output:758d0d40] %[output:410de2f6] %[output:6584ccc7] %[output:5c29f7e3] %[output:09333621] %[output:18f59dd1] %[output:428b866d] %[output:655bd528] %[output:733b13ca] %[output:26b05ffb] %[output:714859bd]
plot3(Y(:,1),Y(:,2),Y(:,3),'color',blue,'linewidth',1.2) %[output:29cc3ebe]
[~,Y] = ode45(@df_pos,[0,40],p0); %[output:47efeee3] %[output:8e819d84] %[output:95d69b61] %[output:3e3e19fb] %[output:127672a9] %[output:99e43c17] %[output:561b6d53] %[output:5dc51cc4] %[output:8a9b98cf] %[output:2a6a1ce3] %[output:674263b5] %[output:7212788e] %[output:8ea23afc] %[output:4f71dd9e] %[output:2eee1378] %[output:451e202f] %[output:5f918b4b] %[output:4cf89415] %[output:3ac229b2] %[output:8f9f400e] %[output:29ee74c8] %[output:2d6ad4ee] %[output:472ea4af] %[output:57567243] %[output:434837c6] %[output:83bf1008] %[output:2ecb41a0] %[output:941b91ce] %[output:03734798] %[output:5ec3dd01] %[output:1c0dae8c] %[output:8d79adaa] %[output:7980d3bc] %[output:8dfa22f8] %[output:2ec06700] %[output:425a6244] %[output:46d9f801] %[output:5ab81d39] %[output:1b0611de] %[output:4c83fc86] %[output:33239548] %[output:0d02d509] %[output:644cb82c] %[output:006a60cf] %[output:4b76eee0] %[output:8e0a8ff8] %[output:5bcf3aac] %[output:9fad7fd6] %[output:032a4ffc] %[output:17289564] %[output:884f0a9b] %[output:4a7a076d] %[output:7d3d09a3] %[output:276cdb2e] %[output:43cf566a] %[output:58438908] %[output:8ce94e40] %[output:47918bf2] %[output:7c645833] %[output:6c96c3c1] %[output:5bcf194e] %[output:7e189efe] %[output:54a84d90] %[output:2d61c8a4] %[output:7c287cec] %[output:0df6fd2a] %[output:6fe26d37] %[output:2fab11c0] %[output:8ed32339] %[output:69883acc] %[output:16624ed0] %[output:64baab91] %[output:9c5b4127] %[output:9c41febf] %[output:1a0571fc] %[output:6043e827] %[output:54d49fe2] %[output:5f203238] %[output:706a2351] %[output:33ba26cf] %[output:0fe54b48] %[output:7454e599] %[output:23a51775] %[output:4a38a6a9] %[output:1bf7d3e5] %[output:24f256e4] %[output:6dc4e3ae] %[output:6793b544] %[output:375ca9b7] %[output:8be66e15] %[output:4a55ca3c] %[output:272b951e] %[output:48ea53f5] %[output:52370ad1] %[output:5ec5508e] %[output:36533eed] %[output:67895ac1] %[output:1d19d2bd] %[output:8714bfc9] %[output:1b328024] %[output:37eee841] %[output:3d7e41a8] %[output:325d7d22] %[output:0c0c93a2] %[output:9391f6c2] %[output:44d9deb5] %[output:79c00472] %[output:42ef4d76] %[output:875f7b81] %[output:4deb7c0f] %[output:8f39f4e5] %[output:4593e3b6] %[output:838ec305] %[output:979106b9] %[output:31a1add3] %[output:200445f7] %[output:1a411f21] %[output:0b3f2b15] %[output:447787cc] %[output:016584d6] %[output:09457dbd] %[output:8a52dceb] %[output:5d11397a] %[output:3382d6da] %[output:11878de3] %[output:271775ea] %[output:0212c08a] %[output:0d96b856] %[output:2bb0a28b] %[output:3377f262] %[output:8a0e4c5b] %[output:59a92f2c] %[output:19195fdd] %[output:6ea07250] %[output:67567710] %[output:8c2a5480] %[output:4edcc433] %[output:5910e228] %[output:7f4d5938]
plot3(Y(:,1),Y(:,2),Y(:,3),'color',blue,'linewidth',1.2) %[output:29cc3ebe]

view(80,20) %[output:29cc3ebe]
axis equal; axis off; %[output:29cc3ebe]
axis(axisrange); %[output:29cc3ebe]
hold off %[output:29cc3ebe]
print('-dpng','-r300','fig2.1.png'); %[output:29cc3ebe]
%%
function vf = df(~,xyz)
x = xyz(1); y = xyz(2); z = xyz(3);
num = [1,1,1,0,0,0]
alpha = num(1); beta = num(2); gamma = num(3);
d1 = num(4); d2 = num(5); d3 = num(6);
vf = [x*(alpha-alpha*x-(alpha+beta+d1)*y+d2*z);
    y*(beta+d1*x-beta*y-(beta+gamma+d3)*z);
    z*(gamma-(gamma+alpha+d2)*x+d3*y-gamma*z)];
end


function vf = df_pos(~,xyz)
x = xyz(1); y = xyz(2); z = xyz(3);
num = - [1,1,1,0,0,0]
alpha = num(1); beta = num(2); gamma = num(3);
d1 = num(4); d2 = num(5); d3 = num(6);
vf = [x*(alpha-alpha*x-(alpha+beta+d1)*y+d2*z);
    y*(beta+d1*x-beta*y-(beta+gamma+d3)*z);
    z*(gamma-(gamma+alpha+d2)*x+d3*y-gamma*z)];
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:8cbedc38]
%   data: {"dataType":"symbolic","outputData":{"name":"vf","value":"\\left(\\begin{array}{c}\nx\\,{\\left(\\alpha -\\alpha \\,x+d_2 \\,z-y\\,{\\left(\\alpha +\\beta +d_1 \\right)}\\right)}\\\\\ny\\,{\\left(\\beta -\\beta \\,y+d_1 \\,x-z\\,{\\left(\\beta +d_3 +\\gamma \\right)}\\right)}\\\\\nz\\,{\\left(\\gamma +d_3 \\,y-\\gamma \\,z-x\\,{\\left(\\alpha +d_2 +\\gamma \\right)}\\right)}\n\\end{array}\\right)"}}
%---
%[output:98e1ddc4]
%   data: {"dataType":"symbolic","outputData":{"name":"F","value":"\\left(\\begin{array}{c}\nx\\\\\ny\\\\\nz\\\\\nx+y+z-1\n\\end{array}\\right)"}}
%---
%[output:5032ddc0]
%   data: {"dataType":"symbolic","outputData":{"name":"K","value":"\\left(\\begin{array}{c}\n{\\left(-\\alpha \\right)}\\,x+{\\left(-\\alpha -\\beta -d_1 \\right)}\\,y+d_2 \\,z+\\alpha \\\\\nd_1 \\,x+{\\left(-\\beta \\right)}\\,y+{\\left(-\\beta -d_3 -\\gamma \\right)}\\,z+\\beta \\\\\n{\\left(-\\alpha -d_2 -\\gamma \\right)}\\,x+d_3 \\,y+{\\left(-\\gamma \\right)}\\,z+\\gamma \\\\\n{\\left(-\\alpha \\right)}\\,x+{\\left(-\\beta \\right)}\\,y+{\\left(-\\gamma \\right)}\\,z\n\\end{array}\\right)"}}
%---
%[output:35b79825]
%   data: {"dataType":"symbolic","outputData":{"name":"sp","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n0 & 0 & 0\\\\\n1 & 0 & 0\\\\\n0 & 1 & 0\\\\\n0 & 0 & 1\\\\\n-\\frac{\\beta }{\\alpha +d_1 } & \\frac{\\alpha }{\\alpha +d_1 } & 0\\\\\n0 & -\\frac{\\gamma }{\\beta +d_3 } & \\frac{\\beta }{\\beta +d_3 }\\\\\n\\frac{\\gamma }{d_2 +\\gamma } & 0 & -\\frac{\\alpha }{d_2 +\\gamma }\\\\\n\\frac{d_3 +\\gamma }{\\sigma_1 } & \\frac{\\alpha +d_2 }{\\sigma_1 } & \\frac{\\beta +d_1 }{\\sigma_1 }\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\alpha +\\beta +d_1 +d_2 +d_3 +\\gamma \n\\end{array}"}}
%---
%[output:12b4514a]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:1a9f4afb]
%   data: {"dataType":"symbolic","outputData":{"name":"vfi","value":"\\left(\\begin{array}{c}\n-x\\,{\\left(x+2\\,y-1\\right)}\\\\\n-y\\,{\\left(y+2\\,z-1\\right)}\\\\\n-z\\,{\\left(2\\,x+z-1\\right)}\n\\end{array}\\right)"}}
%---
%[output:8cda668a]
%   data: {"dataType":"symbolic","outputData":{"name":"spi","value":"\\left(\\begin{array}{ccc}\n0 & 0 & 0\\\\\n1 & 0 & 0\\\\\n0 & 1 & 0\\\\\n0 & 0 & 1\\\\\n-1 & 1 & 0\\\\\n1 & 0 & -1\\\\\n0 & -1 & 1\\\\\n\\frac{1}{3} & \\frac{1}{3} & \\frac{1}{3}\n\\end{array}\\right)"}}
%---
%[output:0341fa2d]
%   data: {"dataType":"symbolic","outputData":{"name":"J","value":"\\left(\\begin{array}{ccc}\n1-2\\,y-2\\,x & -2\\,x & 0\\\\\n0 & 1-2\\,z-2\\,y & -2\\,y\\\\\n-2\\,z & 0 & 1-2\\,z-2\\,x\n\\end{array}\\right)"}}
%---
%[output:1d4f8e38]
%   data: {"dataType":"symbolic","outputData":{"name":"spi","value":"\\left(\\begin{array}{ccc}\n0 & 0 & 0\\\\\n1 & 0 & 0\\\\\n0 & 1 & 0\\\\\n0 & 0 & 1\\\\\n\\frac{1}{3} & \\frac{1}{3} & \\frac{1}{3}\n\\end{array}\\right)"}}
%---
%[output:52d4e928]
%   data: {"dataType":"text","outputData":{"text":"------------------------------------------\n","truncated":false}}
%---
%[output:2b7b0eed]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:63163dba]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"lam","rows":3,"type":"double","value":[["1"],["1"],["1"]]}}
%---
%[output:60de2f77]
%   data: {"dataType":"text","outputData":{"text":"------------------------------------------\n","truncated":false}}
%---
%[output:2166523e]
%   data: {"dataType":"text","outputData":{"text":"------------------------------------------\n","truncated":false}}
%---
%[output:17d781c8]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"lam","rows":3,"type":"double","value":[["-1"],["-1"],["1"]]}}
%---
%[output:8f519491]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n1 & 0 & 0\n\\end{array}\\right)"}}
%---
%[output:0666ac36]
%   data: {"dataType":"text","outputData":{"text":"------------------------------------------\n","truncated":false}}
%---
%[output:015dfcca]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"lam","rows":3,"type":"double","value":[["-1"],["-1"],["1"]]}}
%---
%[output:49bd0171]
%   data: {"dataType":"text","outputData":{"text":"------------------------------------------\n","truncated":false}}
%---
%[output:89d2ef0b]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & 1 & 0\n\\end{array}\\right)"}}
%---
%[output:81a7e044]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"lam","rows":3,"type":"double","value":[["-1"],["-1"],["1"]]}}
%---
%[output:131321f4]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n0 & 0 & 1\n\\end{array}\\right)"}}
%---
%[output:41a89f9f]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{ccc}\n\\frac{1}{3} & \\frac{1}{3} & \\frac{1}{3}\n\\end{array}\\right)"}}
%---
%[output:65ad05ee]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"lam","rows":3,"type":"complex","value":[["-1.0000 + 0.0000i"],["0.0000 - 0.5774i"],["0.0000 + 0.5774i"]]}}
%---
%[output:678614a5]
%   data: {"dataType":"symbolic","outputData":{"name":"vfi","value":"\\left(\\begin{array}{c}\n-x\\,{\\left(x+2\\,y-1\\right)}\\\\\n-y\\,{\\left(y+2\\,z-1\\right)}\\\\\n-z\\,{\\left(2\\,x+z-1\\right)}\n\\end{array}\\right)"}}
%---
%[output:20369ba9]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"spi","rows":5,"type":"double","value":[["0","0","0"],["1.0000","0","0"],["0","1.0000","0"],["0","0","1.0000"],["0.3333","0.3333","0.3333"]]}}
%---
%[output:48a428ce]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:03c5b764]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:0eb999ef]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:4ca30baa]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:92245277]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:0bb6ff28]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:4789eab6]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:3ee060de]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:2f3f835c]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:1bd5d37c]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:29d8bd80]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:41994cd9]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:14e12ba0]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:4ef2c694]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:229b0ec2]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:3533ea85]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:3024aeae]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:778dab18]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:44077167]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:6c8c13ea]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:98431882]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:444e1e57]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:6a0a3cb0]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:695dbe4a]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:7e75c12c]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:3736a73d]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:8da9db11]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:5ce2ced0]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:29ab2444]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:6537de23]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:12c28644]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:1e32ade5]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:431d25eb]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:5ec6c8a5]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:5a243eb0]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:97d69e80]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:4b1b874a]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:3b4fcfde]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:53f1ecfa]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:3baaa464]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:9a9bac92]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:3229033f]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:2b4363d4]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"num","rows":1,"type":"double","value":[["1","1","1","0","0","0"]]}}
%---
%[output:8fac1365]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"szdsukreeqdehtj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0e55d2bd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"scuvdfvujelroif","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9e68206e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ykfzzrzilgtzvni","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:213f36bc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bmsvzsqhjjaidcb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8aa8bcfd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dtbxqaywkzrfgts","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7438ed62]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fushydelrutmjkr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:72dc17bf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ftddgyltkftwydq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:69b81093]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vhyzrgtuxnvwtje","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:228d8c36]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vxhskolbbwdrrcl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0032b4f3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rbyhmhysstumpdv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0c17e45c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jbissbrvecpkmkg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:25cba3f3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zkycrmuhpegferj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3fd518d5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"puerewysddfeqcv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5cb7dbe1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xnsktfiqedesfel","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7ab91dce]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ootcvxkxjdgwrbk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9c3e44a8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"forbevaoxvoyyvw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3b9d0293]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ylpbqjvaahwfotl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2482b4f8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ikfedwawhgvwmtg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0b5e2b34]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xbchotzcacukoxa","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:23ce012b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"brujjzsqwktngtz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7d477538]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"acimoucxqngcqhw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:41a15f57]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qwbhjcovupiepkw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5849ee3d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"aeukpxulehbyxob","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8c00e294]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"egaamayateabfob","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:12980bd2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tosoxvtkwnepvgs","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:920bc24b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sgpglbucwakwkos","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:47e00a74]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vxjydpxigvtmdtl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:76b0d150]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zqehzlmqhezecee","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:298ab721]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mhlnbtvpmnvqzfi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4671bddc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zbxkmotmxtavcug","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:05438899]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vnldwjxqcjpxopp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9daff058]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"prbptcwzzkxnuyf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1819c86d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qcfygeuxrgztnnd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:69cd2354]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zqzxxoiumbxzrdg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7627c088]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iqzkqkzpnitxsfx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:977f7072]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mvywmecymfhaahw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:26f20cde]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hcrhbfvnajznoed","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6d106ede]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"giljucbhvqlnnni","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2a36a45f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"laudncojedpywju","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0dd4615f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iimgvnwnlfbyhxv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1a197072]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hthxflyilpqxvyt","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:99ce0af3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oxjvnxbjzyjxuku","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3ad2a0b0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vgfpgkggevhoigq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2652b34f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hnbfopnvocjydgd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1b32b151]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"usfuatohwnieuaa","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8a256509]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lskthxhdmxuubpf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5e31f127]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tgqypdhscbjexbh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:187fc282]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ochscuezwcksifd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:328340af]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oenekkcmudjqqbs","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8a6c360c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nbyemwwthpzrxau","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:32f384ef]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rjuhcinixwgcawm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:79f830c2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ukxrgdehbfokfep","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0092d02a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"spigxdfxoeyhgpn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2c043606]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tbrbksakwemkjdq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8bf0e7c2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jvsdtechsuqyzzu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2cc86898]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"omuljpwpuvahomb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6025bd09]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fvibmcudjnpkiet","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:90a0d0c5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mlzbzlyfnyskqxw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:499d8f8b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ayeqvatsvqhlmpd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:04674943]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"crvxrgudufiubvj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:79786ad2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jdqmgyjftgdeuyl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:593b8dab]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rrnnnqshzakltie","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4dfb785b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sczvivcyolirhmj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0ea445d0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vteigpflwufsbeu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:06f2268d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qabpkixwnemitkh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:20bae2b7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vjgtynpumyqiepj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8bffd7b6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ugoumcvivprngxf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:44779c46]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"juodgzygwukwdow","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:175c4043]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"njlwkoucljxorwt","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1744e439]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lybkhcpxkwwwhpf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:09310e8b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"okgyyoxifpbssmo","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:172a0f6e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ubugpoonqygreih","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:542fd33f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xgdxnwftiqkrlqk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7f8271f7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lcndtfwlqgxtwot","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:61558be3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nwxhtutmxfsjjah","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:443becdd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cljxoqfomxpxbgh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2ffaffeb]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cwxnehntnftzpvx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:31834072]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cmtmjepkutlsjjf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8da1575f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cfkhccukdyyzskc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9b7e4d16]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ynalmmpuhpaizai","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:08c9c99e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ezhssqcgdtcnswj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:91d75aa8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uyzreoedvhhlnzh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6d70e2d1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"giirdcykyhecjtu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0837221b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fzqjbbqydkefeop","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:26b4ebba]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"szsesvgrytevlhc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9fdb5700]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ccleotfqlrgeqmu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9bed7a84]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ykhmuzlprpnemid","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:545ba9ce]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gqwmzdsfmjpurmp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7d721b90]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fabnkwadsrwyitb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3dc0f81c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fjkrgmvgvgorjsd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:57cabdfe]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cjbgxoobzsdshrs","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2809412a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ncfykbbfvravztx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6d1654a7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oscmnjfpccabrpc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:937e4eee]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"thpbuidvzojfwsg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:894e4de0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ydpqlphkvzciejy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:06076702]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"soxafksadvgzggd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:193ace1a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xpiwoshyfkhwasg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:752c27a0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ktavtkoxwuuclyn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0bf87d0d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cxvcatswalxrsxc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2ce7cd77]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gknwfthhwoghnfy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:349ecc9e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cxvmeqmrflpcudh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6d837fed]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mmdjqxsdczhwxvp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7e03d500]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tjtyudndfawulyy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0d70d53c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zepmozcrdhipxsx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0ff8600d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dueduqdlfygewoa","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3cf68df6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gxaekcilocowvxo","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:60c0e9e1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mdylygqktmbnval","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:90a686fa]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mnjscuexybtfuwv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:38d6e69d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rdehenaidiysrny","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:43df9d4b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yaqmfjemymsprad","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:60eeefc8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"choiduxjqzjrsvz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:52c04a16]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"blkysjybnumwfsu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:58d167ec]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kulckaizsaldycx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8fcd971c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"puoouiopszwncba","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5c4a75e2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"grfuycyzacqpghy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:80d020c8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dkijielcafflalc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:86766d06]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zspyquttogfjfeh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8ed249d5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"avpvyqonvdupjdq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:16d5abcd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vytxzuywjnpvfcx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:254a2ed6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rwynckwtkdqqgvn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:162d7526]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"udctiihakludnbl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:052c4857]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gkzsydoqmmhfrsl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:76bd546d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sjtepabinhhdocv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6e855607]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zoksgyvkcjuqzmw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:75a5e231]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lndbmkfextqyyyf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:69f8667d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lgxikeiqsuoznfl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:91dba932]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xbthjxsjfyyprnz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0f1f0578]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yywgcgcfczfikue","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:728ee629]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cvximvtfftsulyg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:02a63255]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ribvcniivvidsmx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:52a7c19a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bunfconrayrujob","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6bdead69]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rkqkmocrnmikzss","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4a086848]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"msjuarzhyvgckyk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:296f398a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rhmufrylylfacfz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:58d08c81]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lwjcktqjmoyqmnr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:528ed5be]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fiiueyfrsvpotqs","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:791f2c75]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"orjpvekztxrvlzh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3cd0f9b7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"baxsbegojuemxwm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3a0afff4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ywtmidnhalqlocx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5173badd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bpntbvnoxjkvnfx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:288b760c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cnsfpttjuzfmixi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4d3db7be]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kkrrlpbyckmilfn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9b884830]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xxjixronbehiwgf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2c9e3764]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vutpabngwuevimd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4c88a3de]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vwvemulyptbpqrc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2a586522]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"djgwyudpkbcqzqv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:22fc056a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"usrrolatxgxrfrr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8420ffbd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bysbapdsyjnzgpg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8746a4aa]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lijoagozbdllfsz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5dff2e80]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"izetopvvszbyprw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:57f57f70]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rcqsbsecjqikscn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4f48c788]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dokdabyyusirlgw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:806adeed]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ospqvovbchyrvmy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:78557b0c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ifjghtmiqtvdrwr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3a79179a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"efyevquamlrflxr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3dcc7bd7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nwsqkjruckyyasy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8f2a7996]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cstyczrbpcyhowd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3a76b020]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mhjydxzagenwbwl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:07e4bb0a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vnbvksrlibeedpx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:409d7227]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hhcmenvrkismffz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:43a981e4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gmisoaczwjgukbj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:967da8f4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qvgfgqaxyiolyee","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:71dd68fc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wgubonxlfvnvvru","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:54dcedb5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lskukoocasswjtq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3729288e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mnbidnctoeprfsa","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3cdfe3ff]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ihpwxqzqhfwwmzn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4c60964d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yciscbfgmblmqmp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:73921f30]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uiobhgvwsqnwsrx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:98eb2766]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"douqfxboqlcktoh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:181d86fd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fxtojzuaiectgms","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5db6f4b3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bposjnluxllfbpo","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:38bdcca9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dsxhwmbejetqofb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5facba9e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vdncjigghlzgdrv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4491e213]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wxlszugvghiblcz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:69e11642]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yujzffpkjcqvysq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:115999ff]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"giyffyegguuvjzh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8a331cf7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nfcpvrijbnjdvqj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:657675d8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qvvwlvkklszoivr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0bcaf2f0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"edhugayxgywfvea","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2800d1f3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pizjljcbuzbbtch","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:22b1eaef]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qslgspndsotwbum","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4e9a7cfc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gzotiqqwlgwuysl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:748a19e1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xffitvfqdtbxmxh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6fee66d3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gzwylimgjczhbrn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1f7e7fd3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ttnflkxnhdqdmic","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:98aa73e6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tctojttxpovzpek","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:50010a4b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nxgqkestdwiookb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:11b4cf06]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pddkrrnmhqxozfy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7eb4fb62]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xbrjsvxryvkbbtq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:30724f9c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"djmxkuzzhtcojtq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9a959871]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"eluiisjbwhrfnmo","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5ec4f0aa]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ibgdwauhgbkowmy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:79e9f0ec]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jcbwrmncbetfrei","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:838b18cc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tuiwyaeyxshlojy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:85273423]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bwzocqwwuqfcfmw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4867c0ce]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ksmuglxcbzmtebj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6b1c2fe6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hbtniznyrdbujwf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4071a987]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iqtlyfusajzkgnl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2e73f8b1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pfdyqkcctpvxrjh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:59e9086a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"otebvyhohijfpmi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5610f497]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xmzjdwetrhqmfqy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:69fb9678]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jpyjelnovuqovpi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5ee42cad]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rxtzvqdmroaesah","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:696dccc7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vnzbxwhsjstntqy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6319a07a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rhzfczdnalyesgw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2007325d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gnmpilsnraapcjx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:03bd9cec]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vftghyizyvhzgpn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1805b338]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ztmeatmggqygeog","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:91b5a663]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dhxckqajjvjejkm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3ee16d6e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ftoojxdremartcs","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:632b57bf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"igfjkapnbrehzhp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5c839f59]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"quqjdizyspgbjpg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9b396e31]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"salhuecmjujmrbg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:658f5540]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"txfdkjoypzunqqu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1a73b7a9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wpncmwukrbztoup","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:157b6084]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gzqjvucvvilmtqp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:873ecacb]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gfrzkdetwzzrzvb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4807c517]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"odnbraxwnpiiwbb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:27bdb534]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fjcbekbxzhuztch","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8d746b0f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sdrsollmgksepqo","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6c0db045]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"isaywxjoiwywkyc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:996777e0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vwckrbuxhcwlhid","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5a4154b0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iqydupfybttkwsw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7415f0e7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kkytjjjoxqpttlj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:862abcbf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xsrnukuvthkezmc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:057a0996]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fbthfsmmtkrytix","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5432e3f9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pjgpesuzxucsctq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:11622cfd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ceykxkrjimoxzpn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6bd84909]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"knpvocrqjmbzkhb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:174f43ed]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nbwrtfvhikjdlug","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:40e848aa]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mdgflitnvkhfwes","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7a6c7b5e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nzwrxjycedzchdw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1be63df6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"apebqhryjxfmkso","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3bc0e1ce]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qyqrszzcudgwjls","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0d01b568]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zzomlhqxbrypxsu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:72acb21e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"eyrlsusghnrqrot","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:77afd7a1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jirhmkskevyjibc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2e617915]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"icmxcpgpmqvcfpb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9887abbc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"abqbbgenepkyabd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:73aa143c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oylqkavucoclhmm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3e1a8a95]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xhxssvbdnwmbkcm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6fd3cf39]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rzgoxonxlysbrrx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:65f764ff]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dhtpkopcrukzlzr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3cff0c2b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zlqqersemzadjwe","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3105e8c6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"eoaisazrxjhcqsn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:63ff485e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"umxeoekiosnlyav","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:67247625]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"obzgvqcxoifgwjx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:320dc3e4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"srwtozsmawbhqti","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0b3c89f4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zaieyeckezbyaol","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4790f2c4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cudirnzcpnzlhzy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7163a3c1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vdmvugxqfqisklx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:297b940c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ysyhqrfisuarcns","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:029a9f03]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pzqciqwrqiyqcin","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:81578a73]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xukxlzzqshujnum","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9068849d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wwahoystujuemst","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:40e226e7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hkeqqfkqafgrjtl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5b20dded]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nzxgpnqrikuzdyp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8769e9ad]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"awpnyylxmprjajk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:80023480]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dvsptxpdvsgvrth","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:15b51bd2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fzlkzqbyuatcyat","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:21c5a378]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vzaucebkwtqolsw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6510b820]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bsgiwsqbggvmozg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:168f3912]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dwrfwfommtedykv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:04e959d6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ansqjznqdkrjkje","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:54fc490b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wrdadosnvpicgkn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0a917969]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dyqokfjtlgublsn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:17692bf1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ztwygurfnduclcv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:57170966]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bgkzusrxzuovtrc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:83157faf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tsuijknvvzihgfa","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:664da028]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"axvarrycwvpfkkc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3b07c81c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wyqnaqmanevcwsr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3d09fde4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tbhudzresulliet","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1eb1ef8e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tjyeqadeczmjcxi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:613ec346]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"aclorjvinpmfrgu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:538ee531]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mhzsuyehysydxaz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:259d9abe]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"npklvbhtmpzhwuj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:096c9d0c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"auimlfkjmuxphpr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:322d5ac4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"saubjscnugthqln","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:62fe6b0a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"psbwqcvbbngnzrl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8a774e59]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"apcozjyfglktijp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3e4f288b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bainlcsewhonbuu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7750f67d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dwiyhbwvkdjcidr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:743cee12]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"twiktjcdmcbcidu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:61fdd38d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xtfwygbhyzsbrue","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:78ccc79d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hsgaaamrdnucuaa","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:031a10f9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lnforflxftwcgof","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:075c14d5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rrwekbugqfdfyhl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:331b52b2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cxmcnbeijfxzypt","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:271daa25]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xdybshgmmxficeh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2345ac56]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"unrzopgswrbpdtj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0ca1c8e9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ibabwzfpcdpkesf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:492118ce]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zhucphencyzkhew","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6e84c9b1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bccgqchdhswvayv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:855260fc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"skvnvdfsheuqqex","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5b585b0c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vtzuzydnbjruaky","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0229dd8c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rbfpyemxzpvuvkr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6c041a82]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pojdckckvebvwbt","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9993e8bb]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"inpynnultvxfhzh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8f7afa11]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"piqsioxyokcrequ","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:86e82cc5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mppwlusuknpnnpi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:86cca90c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wzhvtnmbbpygsjb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:58d3010e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"emsfdhdstghjdqr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2faab95e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rmmxvgxmybmpwls","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2e645e16]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vmcfjybfbtepbhp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:68c0d57e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mbtmdhliqhlabxb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:70a18445]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zvddrszafsiixts","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:178883dd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"erkxkznmfzadlst","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9831720f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vcqgyrnnykztotj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:54e4a787]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uzhukczimapvojb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:996ba7af]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sfcuckbcxdflihd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5e41a2e7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cszdcilstzlnyqk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:00ce9e37]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gbghqeiwrlsyijl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0b8e9600]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kpqtdytsafpbstr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:840fa490]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sggxqrbvjrdhhvs","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9773f19e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ieckimsftfahavj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8d613736]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ayhvsqcmyopducw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:81572fc1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"apksvhpmajisbyn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:56784ed5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nhilhbdagawvnzi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:575b9f6f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tsagqhmphokgbgg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9bbfc1e0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jxfygvzowyruqpr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:573c6ea0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wcrdgvpzktsetcy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1c7d4558]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rfekscvjvkuwhdi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0930483b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kysnshmvcplhuha","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1157a434]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yzcalrcrjvvgyhf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3c18aeab]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ixgrfpudkryzotw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5c5814ee]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jbbrxnhtnqebkif","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7964d2fe]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"trobmeohrcfgjlf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:97b9b513]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"viafzjxjqoyislp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:111287c5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ovzinicjsvmyntx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8d7d6ed9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"spmlczgqxluxpde","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:860c331b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zzkkuazpxrstljb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:25f2f9cc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wqwzuximvvjblxr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:054c0ea1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pmucxcaorwuibea","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2e87fcbb]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xkjwmovhwxeuwgi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:489ab95c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iiizybmblupigca","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4a595c2b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"icraumhcdupurof","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:964f0080]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zrhgyvymetsbckl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2eae40fc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"omrebglseezdzcn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7292aaa6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"spxougglahtqzht","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2d9effa8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ejgoupehlnrjtrt","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:32c37d9f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"knugaezfuqqdgsj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:27a486b0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jvonhpkrtgkimnm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:77ccb615]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hdcwcssrcaehxsx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:71186c69]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iwvbwirdfgkwofk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2b0b8a30]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nqqrasfdzzwsmpl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:32086bb9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zrepgfvmirgdabt","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:29277dbd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"liqkljwtirmeuvv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3502ea4a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"unrlzggqwgdufhj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:74b7d7ac]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cdfjvieltcbehfe","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:61adffe2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wvxalypspfirenv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5f3fdf94]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"doiqkdfaqvctquu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7f42369c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"whlvcnxslgesxwv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:05f9d35a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"neiirwruastdeqk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:38f90cc0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rxnhixhatggixdj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:721db9dd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ppduqjrzgdvrxre","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3a965c1f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"whjggthiclezzco","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6325deb8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ebyewftgluowxqv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:067775cd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qtdigmuiprrxrgm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:103f627b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sunahzhigylsykr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6c5d2ae3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jesoxugeynimuws","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:85013286]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pivepauztvgjdjq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9f509151]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ieymeyzggsqxyot","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8bcbb9f7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zjyttgvugqcgoor","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3e1dc70d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pzbwkhgwkxanikf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:396d5105]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"djhobzjewuzfsgu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:31c8de88]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"csmbrsbvjnhqfdc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:75d3f1c7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dvvmybxxhdmqiax","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0d5aa878]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ilzvxemakgnavbv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1996b230]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nbzbltaghjeywie","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:82e06869]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"afykdtmzbzhxxkb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:62fdfd67]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hyexprkwlwyjofm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1c3a9219]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mcpdguzbeqbwyve","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9f5269c2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jriggprinuzwnim","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:19803bc1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"huiokvcntilzwft","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2563ab05]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ybtuzvusrrzqffg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4ab734e4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uiuryhrtlytoddb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:86d8f36b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vrvlksebzezwanv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:830bc724]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pggrdidhrtvbfmj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:94bc8491]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xlcjilwdpbhylwb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:99cbf1d5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bmnwibbcynaeljg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:03316099]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nghckqnuvpjnxbp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:26cb6bfe]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yyvlyksxrsmkdcc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9896a64d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bmhluurmbdugkji","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:886524af]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ymyasmsubgfbdam","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3db67aab]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"escerkdvgqgheiu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:85e54546]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qyskbxsxzlzbqpt","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:337c6343]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hxiasewhdycholc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0a549124]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"htdlxvylhewuszz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:24969e3c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hokareiwerwpwzu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:366fdbc5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ahjimlihpqbkxre","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:96ba7690]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fekntkkhzibxfzy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:32688f23]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fwkghewycbgvupk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7a174383]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"foltwqbnkkecbuz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3a3b4447]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oxaxcrmbtdyvmtk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6f060927]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"olmhbpfyncyrjok","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:12f89498]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ckacjaonokpccus","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:717d6dc1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gdwevjxxcrvkwyt","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7adede2f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"euwjxwlkkttfdcc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:98319ac6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gwdwcxqxmublsva","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:03f75bb0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uxbzkvrbtlstxxz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9a4527e4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pxlsvhpdhtoulsd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7abe1f36]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gmaacnbrmzjizax","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:958d5cf2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iqdkadxuetzdmzu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9d16d211]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kivxlzemsjppccl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9fed4249]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"leoivizxkhibqvq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:393dab2d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kvljgcdqcdsytns","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:69aafaf9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ywlmlllorxfgdsw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0dfb7868]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"udhrctbdtklpfux","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:56aaa37b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"btbwklwxovqhigr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:59d5ca29]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lmgicuhpayyhjhs","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3e97accf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hyrhdnagxhiqmml","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2d4449c5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"blcadarnsrghsvw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:368e07aa]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"axqpjfeyrkcapzx","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6feaeb46]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"skpgbnokuwnxsff","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:85ba7983]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ckvpwzvyzsahdpl","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:49cd605d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yxclzlqmuypsbrt","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0a4d4382]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zorfopphlqlcubq","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:008394c8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"knswmripzlyohcf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4b8f7a17]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kaqlltzbdrgtiyl","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:630aa356]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hgkcimmfggefnpp","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:323fb556]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uavdqthtyuziutz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6855e9ad]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"aohuemqcatrmyfk","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6c59756f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vleablbfubajzns","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:69570c30]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vuahokvaludokih","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7f2dd004]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xzqxcanimzoinbt","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:60edf16e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"psjfskqotyuagaq","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7b7980cd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"exzctwpsasxyogc","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:9fa91884]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dkhnwniipkyrjmy","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5fe36a5c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nwuhvzbkcikzifa","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2357fd62]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ihdpwxduofzsvgl","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8bdefd2b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sqqjspacxzrbxbi","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5d5b39ff]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"levtgnuiumbyppl","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5fdd333c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hfvdxryxagetylv","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:67cfe10b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ngfuhwsdkijjayr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:961911c1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ajggjkalhuoyikt","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5fbcfd8a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uphzyquputmblcn","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4ee6508d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pcbrmfxruwvcqiz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:71cfed9d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"etweffdqagqioca","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6288a73a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hksyjpeshcuimtk","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:9ba075ef]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pfidvmksnpqerit","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3be775b8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gipeildoaeiodst","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:47d0543c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"crwublthbwasdaf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0352f9f2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ntlmjnirscnqour","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:98c62e4f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rpwepvyscjsdufc","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:35c6f0e9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jvwqkgjynlnefzk","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3e42848f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nggqhzgnfdvlcct","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:17bbf2c5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"epvgzjefendlvsz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:890b6cac]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"chfenlvuvxpsbzz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:68bb5fa3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qfvahxxqxbnxeld","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:73a1a15d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wtuttikjemjqous","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4f3edcf2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xyfmjzsagstxqlt","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:31726c99]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pxnwryodqjlybnf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:63bc83fd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ifbkikgkzpkbnux","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6b62f62e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cfwlnehqfcjoxrh","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:468bc253]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mtcofamvdrkszfa","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:996af791]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pyilikginzbwrto","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:01b4315b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vvxjeyzduyxdtcv","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:87524361]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jxnqyrpogojeswr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7dc7fe0e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qqkazyvqbgvudqt","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:758368d4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vqfkisuqdryvtak","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2a336769]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kusmsendzcdwata","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5ae0ff06]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wtzpdvzrekbisuf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0f572356]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cvwfshrpfwkyjtw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7b038bf6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"omkwccqpujxgdzb","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1a0ab645]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"txohdkczfkglbku","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:603ae915]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kwzbtlmrurghdhx","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3fe04ca0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cctjgettwjtaomk","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:751bd134]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mvhotbihcpiwifu","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8cedb0b8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zrollgvfbaujoal","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1ca336f1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"splvvnrnurxizdk","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3a52a283]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iqcetdmlzrvknpb","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:447875c5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rfgbbrjgcnexrom","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:70223ab6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zjrchxzmbeaapfe","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:349c4c29]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zzpzneccxttdsyy","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:58f2c191]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tvoehlfrbplbost","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:48100f33]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dgxjgndzigjtgsf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:560398d4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bvokocwgkpibrsl","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:676dcdce]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hrkwlbxfppqrvqb","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:62549ce2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wiionaqkcercvqg","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5789374a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"povjqrlulaztkib","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:29a86abe]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"emlogjsuhydfgyj","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:014f258f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"plrkqfsxsiolacl","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:537dba83]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cjjsopdigyghbeb","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:687e1cac]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mllligmjuudjjuv","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:605e4f11]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"aaflkxuaypgrmob","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:82fcb9a2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fbepyzzihaihyce","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:164a1b9b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iqkjeethokgichh","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8680788e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dbpvazmzefyvttd","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:9d007e56]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"konnzzrwkorufkw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1411562e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hnmkobysbwqqscy","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4c9a98db]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wektfqccnrxfrgf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:99053e30]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cpkhqzlkvimvhnj","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:94b32cc8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jhixsfzwionmkbw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2e2989b8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sxvziynyusiaoyd","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6d7360bb]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mwyyqogppxqdqmg","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3e163704]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hpgaydwzniwhevm","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6299b67e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ywmlhrsxpfgsutf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:30dce336]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wohlxmlzedtluaw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3977cbf6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tbrotqoigwexpgw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:91984ce2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qlhjsngitrppauo","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:85e9b295]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nfvgbctbyelcxcv","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8395638a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"avrbfhuwixpckzp","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4b7c1319]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kchrljqiousdzor","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:38912115]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xsgyeurtkdviuee","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8d275658]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vutgyzlqmjpyckr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:04dc5795]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ofjminshhpwngzk","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:28f8d4eb]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pdjxujpypkmabhq","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:95761ab5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ginhuogwahegjty","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2f32a4be]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rutgdhdaurklxyw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0e01745a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jeagowqblkztosi","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:52c7210e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ixxanzkuoiwmhhc","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0adfd303]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vjujfvxzbhcqbuf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8976b4ba]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vwlgvwtpjonjypu","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:20e85172]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lzvmtsmxkvlddim","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4c42c122]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yqwbxczlfwiskeb","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1c2ceb15]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yeisldrbjugmafw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5c8b0d75]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cnxssgjyeuizyoo","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6efe9c59]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wamnpmshhxbvxnn","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:215dd0f7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sfobfzdjsvbcbku","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1500daaa]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"npitzscwzpuhmdh","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:774cbab4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qcbrcjpmyoddlzy","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1acfea74]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"msbwfntkvnikuwm","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3742a461]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"madjxdckcxcpqwi","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:22d724dd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dpsjaaowsnckeyb","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0c96a232]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cbnrmvgebebivrx","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0f9564f2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"efzbylbzewaggis","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:64cc28c1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"igukzjfgvymycpg","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:80d7af2f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ajpkhlxtzyeyirj","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6acc7d40]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"edshzdvxoisavlr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:99178dba]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nxaopfpkdpdhrlu","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:9baf54c6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"axjugqrunfutlrp","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6cbe9aa7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gwdkrmquwzkugbd","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4fec3391]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"drvekyaqtqywlux","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:370ea9d6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"emwcawmrzhwborb","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:96c6f664]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zbesdaosoerzwis","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:36a74d55]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bffcnqpupkqcpne","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:63540d10]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xageopkhjnhtopl","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4a3572f7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"drcsclzolejccdo","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:31ecf745]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rkdnmasjsazhbbp","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:17a1ca9d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"coxwjjfkrzubuml","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:9dc2287a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"brhxyrugjdrpbyf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:026938f3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cgnjgtdfmgmzvto","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:23bb3469]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zpxhbttdaxkxnmx","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:17b7f131]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vulustpzwrakeyy","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:66a36a1f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xnmfjbfgzyyrfwm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:780f318d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ibdrnjgsgxsmyzh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2e891c64]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ohoaasriyzoqcyv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:73b549b2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"banujlclpiicrri","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2ccb067b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zlgxawangwqbees","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9a259f81]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ibevorlwllkxixg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:72f55d9b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uowunsmmhvmfkbn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:89e24c96]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lvklmoirxgdktgr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3cb5ab61]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uahjsnvflibyysw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:87457d59]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zsyjpxdjfznupaa","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:49fdd6af]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mmobonotjnrbbvr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:603e8026]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ajzouxfnpqtqezc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6750cbbc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ydzzleimtphbtgo","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9c2e7b92]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zlagdtjvmamtmsn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:155dba2e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"philjlsqecwbkrp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:91509444]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wwvugfxiiqkluyg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:19b84c56]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tfvatbagsdvkbhh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:95ced293]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kplsywpilpjczfe","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9163499f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"eicxnwkgdtkicdl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2eb91684]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"srojttrcpktrztp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:57edf667]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ttadoltubfqrsvt","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0d7079cf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iqefijkzglloioq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:10cda76a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hmummcxaedaxmgk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3bd85b4a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kzrksuyzjwqvdwb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:37cd5761]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lhghtxcqnxvicuc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4f56fa4d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qmugwrksckuippb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:49de5903]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ojaagxcmbfsawnk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8ab9444d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ajbbtsxnbozkdkh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:51c5d869]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bfmkoslvxxhwnea","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:14038319]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ebapkrrlomlqpnu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3b5bd3f1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"huiuejvxdgcdznf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:297df7ba]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"eixqxafuojfwzvs","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:615b3336]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cmeojdjlmhmxady","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1bb6b617]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mqnbcbyiaebcjea","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3522fa7b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mostothxupvfmlu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:27033efb]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"scbnaeynsnuunha","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8fc6eb8b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xmqwnctoxtvazpl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:91b0ffca]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vryxroccdaehvbr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:97f567d4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zgvdodsteqyiglv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:55e615f0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pfdofiemksqwwfg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8c45b033]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nevheolbljkrzoy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:93ce9353]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ldtvzfzcrdrscyq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:29a973bf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rxtizoxrhyonvnf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0eaf5d2d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nvyycgfesniyydk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:56a241f9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uwrcnorslbxizgw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:40463efd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jjciaikucnczmld","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4962730c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sdvwwvevskhtpkl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9ca36a46]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gmusypuxzsgwrfu","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1930a20f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rmdyfwmtkznvijp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:77dbc6ec]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gbjqejrysdisvee","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:4e78437a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kofosszrizntjag","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6255467b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"svqipfoyzuxgpbd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9d7a9c01]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vwvdaxmlyvnbncg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:239da291]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yefqtyzkgmcelcy","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:835a22e7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zfgvjfxqdvgglsv","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3feb62ff]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ufbtlnbonnsxdny","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8d7c810b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"aaftutygbazndkm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:146b3e42]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ztpstuaulqvuwyd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:7a8346e7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cevwgkanfxxrfqw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:05a0b632]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wciyfnmifnzeqxa","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:56f865bc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"owjpblhkjgenbin","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8f5b2dac]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fqjjghunpvgzydr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:68416f78]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gjpoxfcplcjyzoj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:151afdef]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sggybqffrnagkyg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:55fd7c1b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jtidylpyetkacbc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:73177822]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yceoworzdwcpltz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:730f8e2e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yovfvmgavvlppfs","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:94e58b25]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ftmhkcqcyaieczg","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5fe03690]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"znkwazgjfoetgpe","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:49c96240]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xowhckocblbheqf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:08fd31bd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yzygwcpttwjxbby","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9b33ae04]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wnmxztykysmfxht","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:06bccb80]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ehuygccnvrlfswl","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:18521f56]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xupkyiffdmiowep","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:00718c4b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bhielypudfylxtf","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:8de07be5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dqmnaqlsjxqdvjm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:949214ed]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rbkhdrgcatowabw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:57b07787]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wtsppkhgqcivowr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:24140659]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qvsjlgjanbbozhe","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:33eaa987]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vtewysmqvybhcep","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:97090215]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tjialcfdhpcekvz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:68fd236b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"srranqhwydckaqo","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:228b0167]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"efqrrwccjdayzga","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2b338af4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wdgihmfujmetcgi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:245d99ad]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"efkvzhlkhqxicjr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:50193d57]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mruidnaeibswojc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:530c57ba]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wawsbpxcszchapc","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:747100ed]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xtqtxxfzgyfhqzq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:77d57786]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oggnrvbvupinped","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:21471de4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"epkdpqxasnmicuw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6a97357f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pofonltaedxxhvq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:287a8c89]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vcppzxlptukarmn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:022a3a31]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mvbydmdkmazcdkk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0e47ca2f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fgvyxcteozkrpeb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:24546a6f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"eehnjiynuigszyb","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2e02dce4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wrgcqnigkpinidm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2b1188e9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"weuxjxassepwcdr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1b4dee9e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vclgxtdblhovqua","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:95f7baee]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xvdycmueeuzdygx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:0a43a5d7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hlucirqyewixavj","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:84c26120]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"phtmtpunopidsvw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5c19dc26]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"aqqnuofslwovylo","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:9c8b6f0a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"izjjycmywettduw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:3d50784f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"teutgrgtczizidr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:1f31e725]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tjwqrxcezakznsk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:433f0e23]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"arkrxlfdogkgivw","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:2c975871]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wrjshfxkbtrdlva","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:433a8922]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"evypzfcovhhywym","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:30b3ae83]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oitrbrxtmzzfirq","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:94f26261]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"veacfrpzkbhdpgk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5cc1fc40]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bcvvuwwudfcdgvi","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:28c31db4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rfswenwoavoracn","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:17fc6145]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xmqekjkanlnlkuz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6aba45e6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wfgdmoxldfpaols","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:87e6c965]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"eysxalcpoayezdh","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:99e14172]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lkknyuwilzylvfr","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:758d0d40]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lwpncqyllnojrry","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:410de2f6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"texnlpjfdomwqkz","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:6584ccc7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pvcuixweblxhket","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:5c29f7e3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qjqlckizdamwash","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:09333621]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sfqvifynglpioed","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:18f59dd1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jamejcoyamuksbd","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:428b866d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fotfsztgafqvihm","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:655bd528]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gkmwirnenrymsnx","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:733b13ca]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oejwyytcjrsucxp","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:26b05ffb]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"aynadnykvyfadct","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:714859bd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qznqxivgxvycqfk","name":"num","rows":1,"type":"double","value":"     1     1     1     0     0     0\n"},"version":0}
%---
%[output:47efeee3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qkwltscxzjtqjha","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8e819d84]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"olhgnyvocchixdz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:95d69b61]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yofyfkyukehskyz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3e3e19fb]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hctgwjyorqwkxsz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:127672a9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dwfaxtaltypaezf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:99e43c17]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oxsbbidujwvvrqn","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:561b6d53]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wuqujnohvutcevt","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5dc51cc4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sjovjjlcilatsko","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8a9b98cf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"gyjyrnqqjxjxyuq","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2a6a1ce3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uyucsldnrfkgouq","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:674263b5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wlqhuiwukqqiwpk","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7212788e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cgqkthvpqphguyr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8ea23afc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wboebxggncemyvp","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4f71dd9e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kwsoulgqmsjmdxy","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2eee1378]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kukfoxiwlmoytxh","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:451e202f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zmkxelvumqepdbg","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5f918b4b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"loxzzbkbsebswav","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4cf89415]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tfwgbyvfdoujsae","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3ac229b2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"inbxzjhndkjqpkb","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8f9f400e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jalxtirbhgnnzmy","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:29ee74c8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nerbxjkyzdhabdd","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2d6ad4ee]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bbarnenekbnvnsw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:472ea4af]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"modngqmlmrjgmab","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:57567243]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nywgkvpntaolxoa","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:434837c6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kzxacvprolkkryv","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:83bf1008]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dcpbrtjrkljzxzz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2ecb41a0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ekdcfpdnblfxetn","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:941b91ce]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rzzhgsprbmvvewy","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:03734798]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xjkqeacovdyxnrc","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5ec3dd01]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bfmpriwsycselfz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1c0dae8c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mlgljvzrntpandc","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8d79adaa]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uecqqsdcferxiqk","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7980d3bc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"eqtvltcmlguwuhj","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8dfa22f8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zttlmunqxwabcea","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2ec06700]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pbigvbttmthyxgd","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:425a6244]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lfrkxspjgtdnwan","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:46d9f801]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vzccuswkvwadvbf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5ab81d39]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"svseeltaxqgzxvz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1b0611de]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kzdxsoswyokhhui","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4c83fc86]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oydxquifmwqvvkf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:33239548]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ofajtlzrdbxzsea","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0d02d509]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jxccsphyxywcpat","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:644cb82c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bkprzuwbejhxxsj","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:006a60cf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"srnfkutsnahjccu","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4b76eee0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"frbqrwujzldpdxi","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8e0a8ff8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ylxfpanvobfnzix","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5bcf3aac]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zionnmhecuzdrod","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:9fad7fd6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rakxvxaqsbfwppu","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:032a4ffc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fuylriadyvmpqxm","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:17289564]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"looghunqtrxvbkl","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:884f0a9b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hyxbnegwuygtozi","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4a7a076d]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qbqwjxvzgkmuedn","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7d3d09a3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rtzzfmopxiiobvd","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:276cdb2e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"brnfignrzcboltn","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:43cf566a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"stdylansnlfcbkc","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:58438908]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"blqbaljjqktkiqo","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8ce94e40]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wgaeuhxurmtpjou","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:47918bf2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dadllehrthmwpay","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7c645833]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bxclryxonaxbjaa","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6c96c3c1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hwodgqzdlqkbivo","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5bcf194e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ykyngehebyuwasp","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7e189efe]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"enyftkviluzbxdu","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:54a84d90]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"srurghgfydlgdfa","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2d61c8a4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tvsxjmpexomrwnw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7c287cec]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"edymqivmyhwfibb","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0df6fd2a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xmvdgsrzzxgoayd","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6fe26d37]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ifoikztuoxyrgti","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2fab11c0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zzbdwrtzzdvuzto","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8ed32339]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wjyeuolyfnpemxx","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:69883acc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dlcdpvghralfcuo","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:16624ed0]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"uojnmwldwpbqzpv","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:64baab91]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xiyrksmgcccffvn","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:9c5b4127]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"cfdgcsvnftcryqj","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:9c41febf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hihzjoyrhqwdtzm","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1a0571fc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"byxlqxfbinrkqtp","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6043e827]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xdxpvqlvhjrzbiw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:54d49fe2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ontrqabuzipkrwx","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5f203238]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"lbexvcpkwscousy","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:706a2351]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nvnsxqjvkjxoqcm","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:33ba26cf]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bqmzjhagbnazoug","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0fe54b48]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vriwdlmiucpdrls","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7454e599]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"eiqkbxgrjbtllte","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:23a51775]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"pekcrhdtjprakga","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4a38a6a9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"evmjaswkqsaejbu","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1bf7d3e5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"iohdpgmcnulqanr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:24f256e4]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xvfilocnreetagr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6dc4e3ae]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sotmzirautytrss","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6793b544]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"aryzcokvqszlmrm","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:375ca9b7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"umqrtsmawqkdcnm","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8be66e15]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yfuzuswytjrnaos","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4a55ca3c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bafpazlaycfuxac","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:272b951e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fhwzsdvdmaygyxw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:48ea53f5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"aglromhsoiblbeq","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:52370ad1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"unpfjtguszdzsen","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5ec5508e]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xomnyonantdxdsr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:36533eed]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ygstnzlsmcetmms","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:67895ac1]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qingswzdmwkayqo","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1d19d2bd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"slpkdqtmdznxdou","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8714bfc9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kmxboljmfxzcslx","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1b328024]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zbyjsuczknoikfh","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:37eee841]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ottqxzugkpmyhob","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3d7e41a8]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vlezpccfusgckuz","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:325d7d22]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fswhmkhyhvekjvx","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0c0c93a2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ielocpjhxtkwzgg","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:9391f6c2]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qfewashobgxzhqq","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:44d9deb5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"zelmvedilzlhmoo","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:79c00472]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hhnxuzuzubjstgp","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:42ef4d76]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"bzfnndovwjjwpnm","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:875f7b81]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"tfjchirzaupxpxd","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4deb7c0f]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hhpdlceisbvbsuj","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8f39f4e5]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rcglxdaviuysvyb","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4593e3b6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dcipaknygbtocpr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:838ec305]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"mklfnljynnlqzid","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:979106b9]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ifznsiblnpgoclh","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:31a1add3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"sofajzbsjodtwlh","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:200445f7]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"dgbipretgnoaicc","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:1a411f21]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ukqpxftjykeljyw","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0b3f2b15]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"hoeoikutnffqpod","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:447787cc]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jptrozgibabrbqf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:016584d6]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"twvrkmnwslpubuf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:09457dbd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xcfnofrwvpptcfa","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8a52dceb]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"vphesuiucergkre","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5d11397a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"djmzwdpegoktaqd","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3382d6da]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"aqumufuemovznsv","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:11878de3]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"nppqxfkwzmbtfnp","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:271775ea]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"kvktbzcnmpkedpi","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0212c08a]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"qdznlqkdhicdbsv","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:0d96b856]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"oxjqmfpzqjyzrfq","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:2bb0a28b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"havpwphruaqpjjo","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:3377f262]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"ymxwcqtiwywdxen","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8a0e4c5b]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"wmwabzmxlovirzj","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:59a92f2c]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"fjkqkztavtulzxn","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:19195fdd]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"yvhjdfxhthuurnq","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:6ea07250]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"rpwkdmhrxytfsjf","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:67567710]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"xxigbcyheydbrky","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:8c2a5480]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jiomuixtxnyijmr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:4edcc433]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"jgmadszfyaprutg","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:5910e228]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"idmuywyglrpwwtr","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:7f4d5938]
%   data: {"dataType":"not_yet_implemented_matrix","outputData":{"columns":6,"header":"1×6","id":"suxjeoetumugpww","name":"num","rows":1,"type":"double","value":"    -1    -1    -1     0     0     0\n"},"version":0}
%---
%[output:29cc3ebe]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAABJIAAALACAYAAADIe1jZAAAAAXNSR0IArs4c6QAAIABJREFUeF7s3Xuwnld9H\/rlhOEyNUgcMDM1DaCAe6HhKNnJxJ4pczJiJ\/SPkEw4UmfOAPGUS1wINAbhPXHiFNdY3eFmBXBkNi7EiTEUEnwynZI0wg3ExG1xpNQ+QE3saCNngs3Y1sUXjIPAW2fW47zb7371Xp779bNnOlR6n2ddPr+1X9A361nPWadPnz4d\/BAgQIAAAQIECBAgQIAAAQIECBBYIHCWIMkaIUCAAAECBAgQIECAAAECBAgQSCMgSEqj5BoCBAgQIECAAAECBAgQIECAAIEgSLIICBAgQIAAAQIECBAgQIAAAQIEUgkIklIxuYgAAQIECBAgQIAAAQIECBAgQECQZA0QIECAAAECBAgQIECAAAECBAikEhAkpWJyEQECBAgQIECAAAECBAgQIECAgCDJGiBAgAABAgQIECBAgAABAgQIEEglIEhKxeQiAgQIECBAgAABAgQIECBAgAABQZI1QIAAAQIECBAgQIAAAQIECBAgkEpAkJSKyUUECBAgQIAAAQIECBAgQIAAAQKCJGuAAAECBAgQIECAAAECBAgQIEAglYAgKRWTiwgQIECAAAECBAgQIECAAAECBARJ1gABAgQIECBAgAABAgQIECBAgEAqAUFSKiYXESBAgAABAgQIECBAgAABAgQICJKsAQIECBAgQIAAAQIECBAgQIAAgVQCgqRUTC4iQIAAAQIECBAgQIAAAQIECBAQJFkDBAgQIECAAAECBAgQIECAAAECqQQESamYXESAAAECBAgQIECAAAECBAgQICBIsgYIECBAgAABAgQIECBAgAABAgRSCQiSUjG5iAABAgQIECBAgAABAgQIECBAQJBkDRAgQIAAAQIECBAgQIAAAQIECKQSECSlYnIRAQIECBAgQIAAAQIECBAgQICAIMkaIECAAAECBAgQIECAAAECBAgQSCUgSErF5CICBAgQIECAAAECBAgQIECAAAFBkjVAgAABAgQIECBAgAABAgQIECCQSkCQlIrJRQQIECBAgAABAgQIECBAgAABAoIka4AAAQIECBAgQIAAAQIECBAgQCCVgCApFZOLCBAgQIAAAQIECBAgQIAAAQIEBEnWAAECBAgQIECAAAECBAgQIECAQCoBQVIqJhcRIECAAAECBAgQIECAAAECBAgIkqwBAgQIECBAgAABAgQIECBAgACBVAKCpFRMLiJAgAABAgQIECBAgAABAgQIEBAkWQMECBAgQIAAAQIECBAgQIAAAQKpBARJqZhcRIAAAQIECBAgQIAAAQIECBAgIEiyBggQIECAAAECBAgQIECAAAECBFIJCJJSMbmIAAECBAgQIECAAAECBAgQIEBAkGQNECBAgAABAgQIECBAgAABAgQIpBIQJKVichEBAgQIECBAgAABAgQIECBAgIAgyRogQIAAAQIECBAgQIAAAQIECBBIJSBISsXkIgIECBAgQIAAAQIECBAgQIAAAUGSNUCAAAECBAgQIECAAAECBAgQIJBKQJCUislFBAgQIECAAAECBAgQIECAAAECgiRrgAABAgQIECBAgAABAgQIECBAIJWAICkVk4sIECBAgAABAgQIECBAgAABAgQESdYAAQIECBAgQIAAAQIECBAgQIBAKgFBUiomFxEgQIAAAQIECBAgQIAAAQIECAiSrAECBAgQIECAAAECBAgQIECAAIFUAoKkVEwuIkCAAAECBAgQIECAAAECBAgQECRZAwQIECBAgAABAgQIECBAgAABAqkEBEmpmFxEgAABAgQIECBAgAABAgQIECAgSLIGCBAgQIAAAQIECBAgQIAAAQIEUgkIklIxuYgAAQIECBAgQIAAAQIECBAgQECQZA0QIECAAAECBAgQIECAAAECBAikEhAkpWJyEQECBAgQIECAAAECBAgQIECAgCDJGiBAgAABAgQIECBAgAABAgQIEEglIEhKxeQiAgQIECBAgAABAgQIECBAgAABQZI1QIAAAQIECBAgQIAAAQIECBAgkEpAkJSKyUUECBAgQIAAAQIECBAgQIAAAQKCJGuAAAECBAgQIECAAAECBAgQIEAglYAgKRWTiwgQIECAAAECBAgQIECAAAECBARJ1gABAgQIECBAgAABAgQIECBAgEAqAUFSKiYXESBAgAABAgQIECBAgAABAgQICJKsAQIECBAgQIAAAQIECBAgQIAAgVQCgqRUTC4iQIAAAQIECBAgQIAAAQIECBAQJFkDBAgQIECAAAECBAgQIECAAAECqQQESamYXESAAAECBAgQIECAAAECBAgQICBIsgYIECBAgAABAgQIECBAgAABAgRSCQiSUjG5iAABAgQIECBAgAABAgQIECBAQJBkDRAgQIAAAQIECBAgQIAAAQIECKQSECSlYnIRAQIECBAgQIAAAQIECBAgQICAIMkaIECAAAECBAgQIECAAAECBAgQSCUgSErF5CICBAgQIECAAAECBAgQIECAAAFBkjVAgAABAgQIECBAgAABAgQIECCQSkCQlIrJRQQIECBAgAABAgQIECBAgAABAoIka4AAAQIECBAgQIAAAQIECBAgQCCVgCApFZOLCBAgQIAAAQIECBAgQIAAAQIEBEnWAAECBAgQIECAAAECBAgQIECAQCoBQVIqJhcRIECAAAECBAgQIECAAAECBAgIkqwBAgQIECBAgAABAgQIECBAgACBVAKCpFRMLiLQHYHTp0+H++67L\/z5n\/95+NznPhduv\/328Hd\/93fhIx\/5SHjFK15xxkQeeOCB5LMbb7wxnH322WH\/\/v3h\/PPP786EjZQAAQIECBAgQIAAAQIEahMQJNVGrSMC1Qt84QtfCFddddVmR0eOHAnf+973kj+\/7GUvCx\/72MfCOeecs\/n517\/+9fC2t70tHD16dPPvfu7nfi685z3vCc94xjOqH3AHejjrrLO2jDIGdX4IECBAgAABAgQIECAwVAFB0lArb969FHjsscfCD\/7gD4anPvWpIQYe119\/fbjiiis25\/qpT30qXHDBBcmf77rrrvCWt7wlvPSlL03+fPDgwfD9738\/XHzxxeGtb31reMpTntJLo6yTEiRlFXM9AQIECBAgQIAAAQJ9FhAk9bm65jZ4ga997Wvh9a9\/fTh+\/Hhi8Ru\/8RvhDW94Q3jkkUfCO9\/5zvCiF70ovOMd7whPf\/rTw7333hvuv\/\/+JFh62tOeNni7EYAgyVIgQIAAAQIECBAgQIDAkwKCJKuBQI8FHn744eTRtVtuuSWZ5ete97okTPr4xz8e7rjjjvDud787PPvZz+6xQPGpCZKKG2qBAAECBAgQIECAAIH+CAiS+lNLMyFwhkB8VO3KK68Mn\/jEJ5LP4mNtb3rTm8IHPvCB8N73vjc5N2nRz8bGRvjv\/\/2\/h\/\/4H\/9juPXWW5Mzl57\/\/OeHX\/iFX0h2N\/U9iBIkLVohPidAgAABAgQIECBAYEgCgqQhVdtcBynwh3\/4h8ljbPEnnnsUz1C69NJLw2tf+9qF5yDFHU0xiIpvdJv285KXvCRce+21ySNyff0RJPW1suZFgAABAgQIECBAgEAeAUFSHjX3EOiQwG233RYuvPDC8Oijjyaj\/smf\/Mnw4Q9\/ODzvec+bO4u4m2l1dTX8zd\/8Tfi3\/\/bfhh\/5kR9Jgqd4vtLv\/u7vhquvvjo5nPvVr351cl1fz1USJHVosRsqAQIECBAgQIAAAQKVCwiSKifWAYFmBb71rW+FN77xjeGv\/uqvkoFcdtllySNpkwHJ5Chvvvnm8Kd\/+qfJYdyTj6+NQqYYKJ177rnhuuuuC+edd16zE62od0FSRbCaJUCAAAECBAgQIECgkwKCpE6WzaAJpBeIO5Hio2x\/9Ed\/lNwUD9x+17veNfextu9+97vhmmuuSc5B2rFjx9TOvvCFLyTnLf2Df\/APwvXXXx9+7Md+LP2gOnSlIKlDxTJUAgQIECBAgAABAgQqFxAkVU6sAwLNC8SDtT\/60Y8mA4kHbn\/kIx8J27ZtmzmweMD2Y489loREs36+\/OUvh9e85jWh7+ckCZKaX79GQIAAAQIECBAgQIBAewQESe2phZEQqERgfX09XHTRReHo0aNJ+2U9ivbpT386\/Pqv\/3pyflL8f\/H8pD7+CJL6WFVzIkCAAAECBAgQIEAgr4AgKa+c+wh0QODb3\/52uPzyy8Nf\/MVfhAcffHDzwO21tbXwyle+MvcMTp48GS6++OLw9Kc\/PfyH\/\/AfwjnnnJO7rbbfKEhqe4WMjwABAgQIECBAgACBOgUESXVq64tAjQKnT58On\/3sZ0MMjX7t134tebTtL\/\/yL5MR7N27N7ztbW\/LNZp4ftJv\/\/ZvhxhS\/cqv\/MoZB3HnarTFNwmSWlwcQyNAgAABAgQIECBAoHYBQVLt5DokUI\/AXXfdFd7ylreEN7\/5zeHVr351uPLKK8MnPvGJpPO4G+kDH\/hAOPvss1MPJgZTd999d\/jN3\/zN8N\/+239LHmWLB3e\/\/e1vD8961rNSt9O1CwVJXauY8RIgQIAAAQIECBAgUKWAIKlKXW0TaEhg9Ejb9u3bwyWXXBKe8YxnhM985jPJzqT4E9\/E9vGPfzy86EUvSv4cQ6LDhw8nB2c\/+9nPPmPU\/+N\/\/I\/woQ99KBw6dOiMz\/75P\/\/nya6n5z\/\/+Q3NttpuBUnV+mqdAAECBAgQIECAAIFuCQiSulUvoyWwUCCGQr\/\/+7+fPNa2f\/\/+8EM\/9EPJPbfddlu48MILN89J+tSnPpW8wW302Qc\/+MHw\/ve\/Pzzvec+b2Ud8k1vclfSf\/\/N\/Dtddd1343ve+l1z7+te\/Pgmp+njgtiBp4ZJzAQECBAgQIECAAAECAxIQJA2o2KbaP4F4XlHcJfTUpz417Ny5MzztaU8Lt9xyS7j00kvDu971rvAzP\/MzYRSEPPDAA+FNb3pT+OpXv5pAjM5Juv\/++8PKykr4pV\/6pfDyl788NVIMpuJup\/g2uB\/+4R9OgqVRaJW6ERcSIECAAAECBAgQIECAQKcEBEmdKpfBEtgqcO2114b3vOc9yV8+5znPCT\/1Uz+VBEm\/\/Mu\/HF7zmtds2SF06tSpJFyKu5XizzOf+czk\/KQ\/+7M\/Cz\/7sz97xvWLrOPOp+uvvz5cccUVyaXjO5wW3etzAgQIECBAgAABAgQIEOimgCCpm3UzagKJwFVXXRUOHDiwRSPuOnrHO96RnIs0+fMnf\/InyZvWvv\/9729+NO\/6Rcx\/\/dd\/nTzWdu+994arr746CaT8ECBAgAABAgQIECBAgEB\/BQRJ\/a2tmQ1A4J577kl2JH3+858PL37xi5O3qMU3tE0LkSJHfBTud37nd5LDsZ\/73OeGN77xjXOvX0R44sSJpI0jR46E3\/qt3wo\/\/dM\/vegWnxMgQIAAAQIECBAgQIBAhwUESR0unqETaFpgFCTF85fiGUnnnXde00PSPwECBAgQIECAAAECBAhUKCBIqhBX0wT6LvC1r30tebQtns105ZVXztwJ1XcH8yNAgAABAgQIECBAgMBQBARJQ6m0eRIoWSAeth0fkfvkJz+ZnNV0\/vnnl9yD5ggQIECAAAECBAgQIECgbQKCpLZVxHgItEDg29\/+dti3b1\/4r\/\/1v4Z\/82\/+TbLraPLcpfh2uEsvvTRccskl4ed\/\/ufDD\/zAD7Rg5IZAgAABAgQIECBAgAABAlUKCJKq1NU2gY4K\/O3f\/m0SHn3jG99IZvDCF74wXHzxxeFf\/It\/EU6dOhU++9nPhhgkjf7urLPO6uhMDZsAAQIECBAgQIAAAQIEsggIkrJouZbAQATiY2v\/+3\/\/7\/DRj340fPnLXw7Hjx8Pz3nOc8JLXvKS5BG2eCbSy172svCUpzxlICKmSYAAAQIECBAgQIAAAQJRQJBkHRAgQIAAAQIECBAgQIAAAQIECKQSECSlYnIRAQIECBAgQIAAAQIECBAgQICAIMkaIECAAAECBAgQIECAAAECBAgQSCUgSErF5CICBAgQIECAAAECBAgQIECAAAFBkjVAgACBOQKTb6SLB5H7IUCAAAECBAgQIECAwFAFBElDrbx5EyCQSkCQlIrJRQQIECBAgAABAgQIDERAkDSQQpsmAQL5BARJ+dzcRYAAAQIECBAgQIBAPwUESf2sq1kRIFCSgCCpJEjNECBAgAABAgQIECDQCwFBUi\/KaBIECFQlIEiqSla7BAgQIECAAAECBAh0UUCQ1MWqGTMBArUJCJJqo9YRAQIECBAgQIAAAQIdEBAkdaBIhkiAQHMCgqTm7PVMgAABAgQIECBAgED7BARJ7auJEREg0CIBQVKLimEoBAgQIECAAAECBAg0LiBIarwEBkCAQJsFBEltro6xESBAgAABAgQIECBQt4AgqW5x\/REg0CkBQVKnymWwBAgQIECAAAECBAhULCBIqhhY8wQIdFtAkNTt+hk9AQIECBAgQIAAAQLlCgiSyvXUGgECPRMQJPWsoKZDgAABAgQIECBAgEAhAUFSIT43EyDQdwFBUt8rbH4ECBAgQIAAAQIECGQRECRl0XItAQKDExAkDa7kJkyAAAECBAgQIECAwBwBQZLlQYAAgXlfkmedteXT06dP8yJAgAABAgQIECBAgMBgBQRJgy29iRMgkEbAjqQ0Sq4hQIAAAQIECBAgQGAoAoKkoVTaPAkQyCUgSMrF5iYCBAgQIECAAAECBHoqIEjqaWFNiwCBcgQESeU4aoUAAQIECBAgQIAAgX4ICJL6UUezIECgIgFBUkWwmiVAgAABAgQIECBAoJMCgqROls2gCRCoS0CQVJe0fggQIECAAAECBAgQ6IKAIKkLVTJGAgQaExAkNUavYwIECBAgQIAAAQIEWiggSGphUQyJAIH2CAiS2lMLIyFAgAABAgQIECBAoHkBQVLzNTACAgRaLCBIanFxDI0AAQIECBAgQIAAgdoFBEm1k+uQAIEuCQiSulQtYyVAgAABAgQIECBAoGoBQVLVwtonQIAAAQIECBAgQIAAAQIECPREQJDUk0KaBgECBAgQIECAAAECBAgQIECgagFBUtXC2idAgAABAgQIECBAgAABAgQI9ERAkNSTQpoGAQIECBAgQIAAAQIECBAgQKBqAUFS1cLaJ0CAAAECBAgQIECAAAECBAj0RECQ1JNCmgYBAgQIECBAgAABAgQIECBAoGoBQVLVwtonQIAAAQIECBAgQIAAAQIECPREQJDUk0KaBgECBAgQIECAAAECBAgQIECgagFBUtXC2idAgAABAgQIECBAgAABAgQI9ERAkNSTQpoGAQIECBAgQIAAAQIECBAgQKBqAUFS1cLaJ0CAAAECBAgQIECAAAECBAj0RECQ1JNCmgYBAgQIECBAgAABAgQIECBAoGoBQVLVwtonQIAAAQIECBAgQIAAAQIECPREQJDUk0KaBgEC1QicddZZWxo+ffp0NR1plQABAgQIECBAgAABAh0QECR1oEiGSIBAcwKCpObs9UyAAAECBAgQIECAQPsEBEntq4kRESDQIgFBUouKYSgECBAgQIAAAQIECDQuIEhqvAQGQIBAmwUESW2ujrERIECAAAECBAgQIFC3gCCpbnH9ESDQKQFBUqfKZbAECBAgQIAAAQIECFQsIEiqGFjzBAh0W0CQ1O36GT0BAgQIECBAgAABAuUKCJLK9dQaAQI9ExAk9aygpkOAAAECBAgQIECAQCEBQVIhPjcTINB3AUFS3ytsfgQIECBAgAABAgQIZBEQJGXRci0BAoMTECQNruQmTIAAAQIECBAgQIDAHAFBkuVBgACBeV+SZ5215dPTp0\/zIkCAAAECBAgQIECAwGAFBEmDLb2JEyCQRsCOpDRKriFAgAABAgQIECBAYCgCgqShVNo8CRDIJSBIysXmJgIECBAgQIAAAQIEeiogSOppYU2LAIFyBARJ5ThqhQABAgQIECBAgACBfggIkvpRR7MgQKAiAUFSRbCaJUCAAAECBAgQIECgkwKCpE6WzaAJEKhLQJBUl7R+CBAgQIAAAQIECBDogoAgqQtVMkYCBGoVmAyPpnXu7W21lkRnBAgQIECAAAECBAi0RECQ1JJCGAYBAu0QSBMijUYqTGpHzYyCAAECBAgQIECAAIH6BARJ9VnriQCBlgtkCZGESS0vpuERIECAAAECBAgQIFCJgCCpElaNEiDQNYE8IZIwqWtVNl4CBAgQIECAAAECBIoKCJKKCrqfAIFeCAiSelFGkyBAgAABAgQIECBAoGIBQVLFwJonQKAbAoKkbtTJKAkQIECAAAECBAgQaFZAkNSsv94JEGiBQJEQaTR8B2+3oJCGQIAAAQIECBAgQIBA5QKCpMqJdUCAQNsF0gRJO3bsCEePHp05FUFS26tsfAQIECBAgAABAgQIlCEgSCpDURsECHRaIE2QtGiCgqRFQj4nQIAAAQIECBAgQKAPAoKkPlTRHAgQKCxQJEwSIhXm1wABAgQIECBAgAABAh0RECR1pFCGSYBAtQKCpGp9tU6AAAECBAgQIECAQD8EBEn9qKNZECBQUECQVBDQ7QQIECBAgAABAgQIDEJAkDSIMpskAQJpBIRJaZRcQ4AAge4LfOKcH00m8YsP3N79yZgBAQIECBCoWUCQVDO47ggQaLdAkTApzsx5Se2ur9ERIEBgFCKNJIRJ1gQBAgQIEMgmIEjK5uVqAgQGIJAmTPr2Qw+Fa6+9NuxbXQ0nTp7coiJMGsAiMUUCBDopMBkiCZM6WUaDJkCAAIGGBQRJDRdA9wQItFtgFCrt2LEjGejRo0dDDJFGPzFMOrC2FtbX14VJ7S6l0REgMHCBaSGS3UgDXxSmT4AAAQK5BARJudjcRIDAkARimDQrSIoO42FSvC6GTfHHzqQhrRJzJUCgrQJ2IbW1MsZFgAABAl0VECR1tXLGTYBAbQKTj7qN70gaH8TOpaWwsbGR\/JUwqbby6IgAAQIzBYRIFgcBAgQIEChfQJBUvqkWCRDoocB4mDQrSIrT3rW8HI4dP74pEAMlO5N6uCBMiQCB1gsIkVpfIgMkQIAAgY4KCJI6WjjDJkCgXoG0QVIc1e49e8JdR44Ik+otkd4IECCwKeA8JIuBAAECBAhUJyBIqs5WywQI9EggS5AUpx3PTdq7srLlbCU7k3q0IEyFAIFWCtiF1MqyGBQBAgQI9ExAkNSzgpoOAQLVCKQ9J2m8d2FSNbXQKgECBKYJCJGsCwIECBAgUI+AIKkeZ70QINBRgckAaTSNeeckLQqT4ud2J3V0QRg2AQKtFBAitbIsBkWAAAECPRUQJPW0sKZFgEA5AkWDpDiKuDNp3+pq2LZ9ezIob3QrpzZaIUCAgADJGiBAgAABAvULCJLqN9cjAQIdEigjSBqFSQfW1sLGxoYwqUP1N1QCBNorIERqb22MjAABAgT6LSBI6nd9zY4AgYICZQVJo2HsXFoSJhWsidsJECAgRLIGCBAgQIBAcwKCpObs9UyAQAcEyg6S4pR3LS+HY8ePb84+PurmzKQOLAZDJECgFQJCpFaUwSAIECBAYMACgqQBF9\/UCRBYLFBFkBR73b1nT7jryBFh0uISuIIAAQKJgADJQiBAgAABAu0QECS1ow5GQYBASwWmBUk7duxIDsxO++a2WVOLh3Bffc01wqSW1t6wCBBoj4AQqT21MBICBAgQICBIsgYIECAwR2DWjqR4S9EgKbYRw6S9KyshhlPxx2NuliMBAgS2CgiRrAgCBAgQINAuAUFSu+phNAQItEyg6iBpVpgU\/965SS1bDIZDgEDtAtNCpF984Pbax6FDAgQIECBA4EkBQZLVQIAAgTkCdQRJozBp3+pq2LZ9ezKauDNJmGRpEiAwVAG7kIZaefMmQIAAgS4ICJK6UCVjJECgMYG6gqRRmHRgbS2sr69vma+dSY2VX8cECDQgIERqAF2XBAgQIEAgg4AgKQOWSwkQGJ5AnUHSSHfn0pIwaXhLzYwJEPBmNmuAAAECBAh0QkCQ1IkyGSQBAk0JNBEkxbkKk5qquH4JEGhKwHlITcnrlwABAgQIZBMQJGXzcjUBAgMTmBckRYoy3tw2i3TX8nI4dPjwlo895jawBWi6BAYg4FG2ARTZFAkQIECgVwKCpF6V02QIEChboMkgKc7l2muvDXtXVoRJZRdWewQItEJAiNSKMhgEAQIECBDIJCBIysTlYgIEhibQdJAkTBraijNfAsMRECINp9ZmSoAAAQL9EhAk9aueZkOAQMkCbQiSZoVJ8e896lZywTVHgEAtAs5DqoVZJwQIECBAoBIBQVIlrBolQKAvAm0JkkZh0r7V1XDi5MktvMKkvqw28yDQfwG7kPpfYzMkQIAAgf4LCJL6X2MzJECggMBkkBRDm\/G\/q\/Kw7WnDjmcmCZMKFNStBAg0JiBEaoxexwQIECBAoFQBQVKpnBojQKBvAm0LkqJvDJMOrK2FjY2NhPvo0aPJf9qZ1LfVZz4E+iMgROpPLc2EAAECBAgIkqwBAgQIzBFoY5A0Gu7OpSVhktVLgEDrBZyH1PoSGSABAgQIEMgkIEjKxOViAgSGJrAoSIoedT\/eNl6DXcvL4djx45t\/FXcn2Zk0tFVqvgTaKWAXUjvrYlQECBAgQKCogCCpqKD7CRDotUDbg6SIv3vPnnDXkSPCpF6vRJMj0C0BIVK36mW0BAgQIEAgi4AgKYuWawkQGJxAF4KkWJR4btLelZWwY8eOpEZ2Jg1uqZowgdYICJFaUwoDIUCAAAEClQgIkiph1SgBAn0WmAyXmny0bdxZmNTnVWduBLohIETqRp2MkgABAgQIFBEQJBXRcy8BAoMViGFS3P0Td\/60JUiKxZgWJsW\/d27SYJeqiROoTcCh2rVR64gAAQIECDQqIEhqlF\/nBAh0VaCtQdIoTNq3uhq2bd+e8MawS5jU1ZVm3ATaL2AXUvtrZIQECBAgQKBMAUFSmZraIkBgMALjj7e1aUdW+iiOAAAgAElEQVTSqABxZ9KBtbWwsbEhTBrMqjRRAvULCJHqN9cjAQIECBBoWkCQ1HQF9E+AQCcFRkFSGx9vGwfdubQkTOrkCjNoAu0XECK1v0ZGSIAAAQIEqhAQJFWhqk0CBHov0NYDt6fBC5N6vxxNkEDtAs5Dqp1chwQIECBAoDUCgqTWlMJACBDokkCXgqToumt5ORw7fnyTOJ6b5ADuLq04YyXQDgG7kNpRB6MgQIAAAQJNCgiSmtTXNwECnRXoWpAUoeO5SVdfc40wqbOrzsAJNCsgRGrWX+8ECBAgQKAtAoKktlTCOAgQ6JRAF4OkUZi0d2UlxLOd4o+dSZ1adgZLoDEBIVJj9DomQIAAAQKtExAkta4kBkSAQFcE2v7mtlmOcWfSZJgUr\/WoW1dWnnESqFfAeUj1euuNAAECBAi0XUCQ1PYKGR8BAq0V6GqQFEFjmLRvdTVs27498Y07k4RJrV1qBkagEQG7kBph1ykBAgQIEGi9gCCp9SUyQAIE2irQ1cfbRp4xTDqwthbW19e3ENuZ1NYVZ1wE6hMQItVnrScCBAgQINA1AUFS1ypmvAQItEag60FShBQmtWY5GQiB1ggIkVpTCgMhQIAAAQKtFBAktbIsBkWAQFsEJsOi8d068bN4aPXosbBvP\/RQW4adeRw7l5bsTMqs5gYC\/RNwHlL\/ampGBAgQIECgbAFBUtmi2iNAoFcCi4KkONlRmNTlICnOY9fycjh0+PCW+nnMrVfL2WQIzBUQIlkgBAgQIECAQBoBQVIaJdcQIDBYgTRB0gin60FSnMfuPXvCwZtuEiYNdsWb+BAFPMo2xKqbMwECBAgQyC8gSMpv504CBAYgMC9IitPv8pvbZpUvnpu0d2VFmDSA9W2KBIRI1gABAgQIECCQVUCQlFXM9QQIDEpgiEFSLLAwaVDL3GQHKiBEGmjhTZsAAQIECBQUECQVBHQ7AQL9FhhqkDQrTIp\/79ykfq95sxuGgPOQhlFnsyRAgAABAlUICJKqUNUmAQK9EcgSJMVJ9+GcpPHixZ1J+1ZXw4mTJ7e8oU6Y1JslbiIDE7ALaWAFN10CBAgQIFCBgCCpAlRNEiDQH4GhB0mxkjFMOrC2FjY2NpLCHj16NPlPYVJ\/1rmZDENAiDSMOpslAQIECBCoWkCQVLWw9gkQ6LSAIOnJ8u1cWhImdXo1G\/yQBYRIQ66+uRMgQIAAgXIFBEnlemqNAIGeCQiSthZUmNSzBW46gxBwHtIgymySBAgQIECgNgFBUm3UOiJAoIsCi4KkOKfxa\/p2RtK0mu1aXg7Hjh\/f\/Cg+6uYxty6ubmMegoAQaQhVNkcCBAgQIFCvgCCpXm+9ESDQMQFB0vSCxXOTrr7mGmFSx9az4Q5HwKNsw6m1mRIgQIAAgboFBEl1i+uPAIFOCQiSZpcrhkl7V1aSt7nFHzuTOrW0DbbHAkKkHhfX1AgQIECAQAsEBEktKIIhECDQXoGsQVKcyRAebxtVbFqYFD\/zqFt717SR9VtAiNTv+podAQIECBBog4AgqQ1VMAYCBForIEhaXJoYJu1bXQ3btm9PLo47k4RJi91cQaBsAechlS2qPQIECBAgQGCagCDJuiBAgMAcAUFSuuURw6QDa2thY2NDmJSOzFUEShOwC6k0Sg0RIECAAAECKQQESSmQXEKAwHAFBEnpay9MSm\/lSgJlCQiRypLUDgECBAgQIJBWQJCUVsp1BAgMUkCQlL3sO5eW7EzKzuYOApkFhEiZydxAgAABAgQIlCAgSCoBURMECPRXIE2QFGc\/ft2QDtueVfldy8vh2PHjmx97o1t\/f0fMrBkB5yE1465XAgQIECBAIARBklVAgACBOQKCpPzLY\/eePeGuI0eESfkJ3UlgqoAQycIgQIAAAQIEmhQQJDWpr28CBFovkCdIipOyK+mJ0sZzk\/aurIQdO3Ykf7YzqfVL3gBbLOBRthYXx9AIECBAgMCABARJAyq2qRIgkF0ga5AUA5MYlgiSnrQWJmVfd+4gMCkgRLImCBAgQIAAgbYICJLaUgnjIECglQKCpHLKMgqTJls7ffp0OR1ohUCPBYRIPS6uqREgQIAAgQ4KCJI6WDRDJkCgfQKTgZMdSWfWKIZJ+1ZXw4mTJ7d8KExq33o2ovYIOA+pPbUwEgIECBAgQOAJAUGSlUCAAIESBARJ6RBjmHRgbS2sr68Lk9KRuWrAAkKkARff1AkQIECAQIsFBEktLo6hESDQLYHxMMmOpPm127m0JEzq1vI22hoFPMpWI7auCBAgQIAAgcwCgqTMZG4gQIDAdAFBUraVsWt5ORw6fHjLTR5zy2bo6v4JCJH6V1MzIkCAAAECfRMQJPWtouZDgEBjAoKk7PS79+wJB2+6SZiUnc4dPRQQIvWwqKZEgAABAgR6KCBI6mFRTYkAgWYEYpC0Y8eOcPTo0WQAHm9LV4dpb3SzMymdnav6I+A8pP7U0kwIECBAgEDfBQRJfa+w+REgUJvAKEiKHcYwSZCUnl6YlN7Klf0TECL1r6ZmRIAAAQIE+iwgSOpzdc2NAIFaBQRJxbinhUmxRbuTirm6u70CHmVrb22MjAABAgQIEJgtIEiyOggQIFCSgCCpOGQMk\/atroZt27cnjY0eExQmFbfVQrsEhEjtqofRECBAgAABAukFBEnprVxJgACBhQLj5yR5tG0h19QLYph0YG0tbGxsCJPyEbqr5QJCpJYXyPAIECBAgACBuQKCJAuEAAECJQp4c1t5mDuXloRJ5XFqqSUCzkNqSSEMgwABAgQIEMgtIEjKTedGAgQInCkgSCp3VQiTyvXUWnMCdiE1Z69nAgQIECBAoFwBQVK5nlojQGDgAuNBUqTweFvxBbFreTkcO358s6F4bpIzk4q7aqE+ASFSfdZ6IkCAAAECBKoXECRVb6wHAgQGJCBIqqbY8dykq6+5RphUDa9WKxQQIlWIq2kCBAgQIECgEQFBUiPsOiVAoCsCk8HQop0wgqTqKhvDpL0rK2HHjh1JJ3YmVWet5XIEhEjlOGqFAAECBAgQaJeAIKld9TAaAgRaJiBIaldBpoVJcYSLAr52zcJohiDgUO0hVNkcCRAgQIDAMAUEScOsu1kTIJBSQJCUEqrGy2KYtG91NWzbvj3pNe5MEibVWABdLRQQIi0kcgEBAgQIECDQYQFBUoeLZ+gECFQvkDVIiiPy5rbq6yJMqt5YD9kFPMqW3cwdBAgQIECAQPcEBEndq5kREyBQo0DeICme4xN3ynhrW3XFimHSgbW1sLGxkXRiZ1J11lpeLCBEWmzkCgIECBAgQKAfAoKkftTRLAgQqEigaJAUhyVMqqg4f9\/szqUlYVK1xFpfICBEskQIECBAgACBIQkIkoZUbXMlQCCzQN4gKXZkV1Jm7tw37FpeDseOH9+83xvdclO6MaOA85AygrmcAAECBAgQ6LyAIKnzJTQBAgSqFBAkValbbtu79+wJdx05Ikwql1VrcwSESJYHAQIECBAgMEQBQdIQq27OBAikFigSJI068Whbau7CF8Zzk\/aurCS7weKPnUmFSTUwRcCjbJYFAQIECBAgMGQBQdKQq2\/uBAgsFBAkLSRq3QXCpNaVpFcDEiL1qpwmQ4AAAQIECOQQECTlQHMLAQLDEcgTJEWd8fvsSKp\/vYzCpMmeT58+Xf9g9NgbASFSb0ppIgQIECBAgEABAUFSATy3EiDQfwFBUndrHMOkfaur4cTJk1smIUzqbk2bHLnzkJrU1zcBAgQIECDQJgFBUpuqYSwECLROQJDUupJkGlAMkw6srYX19XVhUiY5F48LCJGsBwIECBAgQIDAkwKCJKuBAAECcwTKCJJi8x5va3aZ7VxaEiY1W4JO9u5Rtk6WzaAJECBAgACBigUESRUDa54AgW4LCJK6Xb\/x0QuT+lPLOmYiRKpDWR8ECBAgQIBAFwUESV2smjETIFCbgCCpNupaOtq1vBwOHT68pS9nJtVC36lOhEidKpfBEiBAgAABAjULCJJqBtcdAQLdEhAkdateaUY77Y1uwqQ0csO4Rog0jDqbJQECBAgQIJBfQJCU386dBAgMQECQ1M8iC5P6Wdeis3KodlFB9xMgQIAAAQJDEBAkDaHK5kiAQG6BvEFS7HD8Xodt5y5BZTeOh0k7duwIR48eTfqyO6ky8lY3LERqdXkMjgABAgQIEGiRgCCpRcUwFAIE2icgSGpfTcocUQyT9q2uhm3bt28GScKkMoXb35ZH2dpfIyMkQIAAAQIE2iUgSGpXPYyGAIGWCZQVJMVp2ZXUsuL+\/XBimHRgbS2sr69vGaCdSe2sV5mjEiKVqaktAgQIECBAYCgCgqShVNo8CRDIJSBIysXWuZuESZ0rWeEBC5EKE2qAAAECBAgQGKiAIGmghTdtAgTSCQiStjqdfcXhM+C+fflPpMPswFU7l5bCxsZGMlJnJnWgYDmH6DyknHBuI0CAAAECBAjEs2BP27tvIRAgQGCmQBlB0ugg564+2jYtPJq3ZLoeLO1aXg7Hjh\/fnGIMlPxXZX++JIRI\/amlmRAgQIAAAQLNCAiSmnHXKwECAxAYhVBdDZKyBkjTShpDpdhO18Kl3Xv2hIM33RRi7eKPMKkfv7BCpH7U0SwIECBAgACBZgUESc36650AgZ4LxDCpS0FS2vBo\/788Nxw4dCycfOzxcOKxx1NXsUuBUjw3ae\/Kypa52ZmUutStutB5SK0qh8EQIECAAAECHRcQJHW8gIZPgEC7BboSJKUJkGJ4NOvnP33twbB+4ru9C5WESe3+\/UozOiFSGiXXECBAgAABAgTSCwiS0lu5kgABApkFuhAkLQqR5gVI00AOHnkkHDn53bB+4lQqr7bvUpoWJsWJ2Z2UqryNXiREapRf5wQIECBAgEBPBQRJPS2saREg0A6Btp+TtPPqr8wMfEYB0kUXTN+JdO2X702FnPYRuDYHSjFM2re6GrZt3775NjdhUqryN3aR85Aao9cxAQIECBAg0HMBQVLPC2x6BAg0KzD51re2vbltcjdSDI9mBUeLJBcFS7HdXR+7I5x47Ptzdyu1NVCKYdKBtbWwsbGRUMQDuIVJi1ZFM58LkZpx1ysBAgQIECAwDAFB0jDqbJYECDQk0LUgqYwQJ02gFMsRQ6VD93xnamXKGEdVJd+5tJSESaMgSZhUlXS+doVI+dzcRYAAAQIECBBIKyBISivlOgIECOQQaHOQ9IL33bblcOwiu5Gm0cwLlEa7nuI18dG3WecptTVQEibl+GWo4RYhUg3IuiBAgAABAgQGLyBIGvwSAECAQNUC42FSWx5tiwHO3oNbzziqKrRJGyhNjme8LlWNrUjtdy0vh2PHj9uZVASxpHsdql0SpGYIECBAgAABAikEBEkpkFxCgACBIgJtC5Km7QKqI6iZFSiNn8k07\/DvOsaYtc7x3KSrr7kmuc2ZSVn1yrleiFSOo1YIECBAgAABAmkFBElppVxHgACBnAJtDJLq2o00SZZ2d9K+L9235bG78XbaFijFMGnvykrYsWOHMCnn70je24RIeeXcR4AAAQIECBDILyBIym\/nTgIECKQSaNs5SZOPtTURzKTZnTTvMO4I38S4ZxV8FCZNfn769OlUa8RF2QWch5TdzB0ECBAgQIAAgTIEBEllKGqDAAECcwRikBR3q8Sf+PhT0+cknX3F4S2jbSqQSRMmdekw7hgm7VtdDSdOntziK0wq\/+tBiFS+qRYJECBAgAABAmkFBElppVxHgMAgBSZ3E+UJBdoUJNV5yHbaBZM2UJp1GHdTQdi0+cUw6cDaWlhfXxcmpV0AGa8TImUEczkBAgQIECBAoGQBQVLJoJojQKBfAmUFSVFldIZOkzuSXvC+27acPdSWECZNmBQNu3AYtzCpmu8A5yEVc41+v\/jA7cUacTcBAgQIECBAIIQgSLIMCBAgMEegzCBp1E2TQVJbHmubRp42TIrXzTqMuy3BWBJ6LS3ZmVTSt4sQqRjkuJ8wqZiluwkQIECAAAFBkjVAgACBuQJlBEmxgza8uW1yN0+bQpfxIkwLlC664Nwz6jTrMO42zWvX8nI4dHjrmVR5Ho8c8q+pEKlY9T0KWMzP3QQIECBAgMCZAnYkWRUECBCYI9CnIKnNu5EmS5Bld9KBQ8fC+olTZ1SxLYHS7j17wsGbbtoyPmFSuq8dIUg6p1lX8Svm524CBAgQIEBguoAgycogQIBAzUFS7K6Jx9u6FCRFo7RhUrx29w13hoPrj0ytZBsCpXhu0t6VFWFShm8bIUgGrCmX8ivm524CBAgQIEBgtoAgyeogQIDAAIKkybe1tSFcSbPwsoRJ8do2704aD5NGB6\/bmTR9FQhB0vx2zL6GXzE\/dxMgQIAAAQLzBQRJVggBAgRqCpJG4UHsru4dSV0NkqJVljApXj9rd1IbwrMYJu1bXQ3btm9PVt3Ro0eT\/xQoPflLKAQp9pXMr5ifuwkQIECAAIHFAoKkxUauIEBgwAJVnJEkSMq3oNIewj0Kn6a92U2YlM++jrscql1cWYhU3FALBAgQIECAwGIBQdJiI1cQIDBggb4ESZO7dNoQqORZVlnCpNj+5JvqRn02Pf+4M+nA2lrY2NhIhjT0nUlCpDy\/DVvvESIVN9QCAQIECBAgkE5AkJTOyVUECAxUoKwgKfKNt1X3o227PnZHOHTPdzar2HSQUmQ5ZQ2TJh\/ra0uYlARdS0uDD5OESEV+G564V4hU3FALBAgQIECAQHoBQVJ6K1cSIDBAgb4ESZM7c7ocJMVlmCdMauujbruWl8Ox48c3f7vi7qShnJkkRCr+pSpEKm6oBQIECBAgQCCbgCApm5erCRAYmIAgqb0FzxomxZm09VG33Xv2hLuOHBlUmCQAKf67xbC4oRYIECBAgACB7AKCpOxm7iBAYEACVQVJkbDOx9v6tiNptATzhEltfdQtnpu0d2UlxLf7xZ8+70wSgBT\/EmVY3FALBAgQIECAQD4BQVI+N3cRIDAQgT4ESTE4GX+sq+uPtU0uvbxhUhsfdRtCmCQAKf7lybC4oRYIECBAgACB\/AKCpPx27iRAYAACgqRuFDlvmHTg0LGwfuLUGZNsMmybFibFAfbh3CQBSPHfJ4bFDbVAgAABAgQIFBMQJBXzczcBAj0X6EuQtPfgvZuVajIkqXK55AmT4ngm32g3GmOTTjFM2re6GrZt354MJz7m1uUwyaHa5ax8IVI5jlohQIAAAQIEigkIkor5uZsAgZ4LlBkkRarx9uo6I2nyTKAmA5Kql0veMKmN5ybFMOnA2lrY2NjodJgkRCpn1QuRynHUCgECBAgQIFBcQJBU3FALBAj0WECQ1L3iFgmT2nhu0s6lpc6GSUKkcn5\/hEjlOGqFAAECBAgQKEdAkFSOo1YIEOipgCCpm4UtEia18dykLoZJQqRyfneESOU4aoUAAQIECBAoT0CQVJ6llggQ6KFAVUFSfMV7PPemjsfbhvRo2\/gSzBsmxTbaeG7SruXlcOz48c0pxvXT1gO4hR\/lfBlyLMdRKwQIECBAgEC5AoKkcj21RoBAzwSqCJJiiBR\/BEnVL5YiYVJbz026+pprWh0mCT\/KWdccy3HUCgECBAgQIFC+gCCpfFMtEiDQI4EqgqTIU+eOpN033BkOrj+yWZU+H7Y9bekVDZPadm5SPIR778rK5hqKc27LziThRzlffhzLcdQKAQIECBAgUI2AIKkaV60SIEBgqsBkMFXHo22Tj2kNLUiKhSgaJk07N6lJx1GYNLnImgyUhB\/lfOlxLMdRKwQIECBAgEB1AoKk6my1TIAAgTMEmgiSdl79lbB+4tTmWJoMQJpcEkXCpDjuaecmNWkZw6R9q6vhxMmTW1ibCJOEH+WsbI7lOGqFAAECBAgQqFZAkFStr9YJECAwN0yqY0eSIOnJEhQNk9p2blLTYZI3s5X3BSdEKs9SSwQIECBAgEC1AoKkan21ToAAAUFSy9ZAGWFSm85NimHSgbW1sL6+vkW66p1JQqTyFrYQqTxLLREgQIAAAQLVCwiSqjfWAwECBLYIxMfbRodtxw+q3pX0gvfdFk489ngyhiYfxWrTMpgMky664NxMw4v3t+3cpJ1LS7WFSUKkTMtl7sVCpPIstUSAAAECBAjUIyBIqsdZLwQIENgUGAVJ8S+OHj0qSGpobRQNk+KwJx8bbDqs27W8HA4dPrxFtOydSYKP8hYsy\/IstUSAAAECBAjUJyBIqs9aTwQIEEgE6g6Szr7iyWDBjqQnF2HRR9xGLU07hLvJQGn3nj3h4E03VRImCT7K+xJjWZ6llggQIECAAIF6BQRJ9XrrjQABArUHSR5tm73oygqT2ngI996VlVLDJMFHeV9eLMuz1BIBAgQIECBQv4AgqX5zPRIgMHCBunckCZLmL7i+h0nj53HlfcxN8FHelxbL8iy1RIAAAQIECDQjIEhqxl2vBAgMXGD8wO2qD9seP8fHo23TF14Z5yXFlmM7bXujW9yZNB4mxXFmCZQEH+V9WbEsz1JLBAgQIECAQHMCgqTm7PVMgMCABWKQNPoRJLVjIZQZJrXpjW7XXntt2Le6GrZt355AxwPe04ZJgo\/y1ibL8iy1RIAAAQIECDQrIEhq1l\/vBAgMVKCpIClyd2FX0v\/74v9rc2X83+tfqmWVlPWI22iwbXqjWwyTDqythY2NjdRhkuCjvGXHsjxLLREgQIAAAQLNCwiSmq+BERAgMECB8SApCXceeqgyhclAo61B0nh4lBaj7JCp7DBp2hvdmvTfubSUKkwSfKRdgYuvY7nYyBUECBAgQIBAtwQESd2ql9ESINATgTqDpMkwo8kgY7J8ecKjWUugrFCp7DBp2hvdmqzBvDBpWugRvX\/xgdt78ptX7zSESPV6640AAQIECBCoR0CQVI+zXggQ6KjAZOCT5ZDieVMeepBUZoA06VxGoFRHmBTH3VSgtGt5ORw7fjyhG52ZdP1zd05dskKkfF9eQqR8bu4iQIAAAQIE2i8gSGp\/jYyQAIEGBfoQJO2+4c5wcP2RTcWmwovRALKESKNQ6ODya8L3Hnw4nHrw4UyroUioVNbh26MBT9uZ1GSYFM9Nuvqaa5LhXfHIs4RImVbW\/IuFSCViaooAAQIECBBonYAgqXUlMSACBNokUFWQFOdY14HbkwFG24OkNOFP1mApTZuT667sXUmx\/TaGSc997w1CpBK\/dIRIJWJqigABAgQIEGilgCCplWUxKAIE2iJQdZC0Y8eO5NGiKg\/bblOQNGs3Up6gZ3yNxGAp\/jx69zfnLp2s\/VQVJu370n3hxGOPbxlrEwHftHpceOz\/S8ZV1mOcbfldrmMcQqQ6lPVBgAABAgQINC0gSGq6AvonQKDVAoKk8sozLbTIGuykGc1fX39j+Mbv3Tg3VMrSb9mPuMU5xDabDpOm1ePyZz7x6ODo3CRhUpoV98Q1QqT0Vq4kQIAAAQIEui0gSOp2\/YyeAIGKBeoIkkb\/cK9qV1IbdiTVFSKNL4cYKN35oetmnqvUhjDpwKFjYf3EqS2ruI6dSdPqcexXXxcOrK2FjY2NzSApDkyYtPhLRoi02MgVBAgQIECAQH8EBEn9qaWZECBQgUDVQVIcctWPtzUdJFX1OFvacpcRKFXxiNto\/Duv\/kqtYdK8UC8ewC1MSruypu9Cind70116Q1cSIECAAAEC3RMQJHWvZkZMgECNAnUESaPp9HVHUhO7kaYtkS\/ufnM4efsdU1dPmt1JVTziVneYlLYWO5eWkp1J8cdjbtO\/cKbtQhIi1fjlrCsCBAgQIECgMQFBUmP0OiZAoAsCgqRiVYqHYE8egJ0mtCnW6\/y7iwRKXQ6T0oZII71dy8vh2PHjm5gxUPKY2xMcQqQqf0O1TYAAAQIECLRdQJDU9goZHwECjQoIkorxT4YXTYdI47OZFnLFz+eNscpH3GLfuz52Rzh0z3e2oJdxZlLWEGk0gN179oS7jhxJ\/mhnkhCp2LeBuwkQIECAAIG+CAiS+lJJ8yBAoBKBKoOkOODx9vv4aFubg6ToP+8Nb7MCpSp3JcUx7b7hznBw\/ZHSwqS8IdJoAPHcpL0rK1vGM9SdSXYiVfI1q1ECBAgQIECgYwKCpI4VzHAJEKhXQJCU33tyx0+bdiNNzmrW425NhUmTB6TH8ebZmVQ0RBImPblShEj5vwvcSYAAAQIECPRLQJDUr3qaDQECJQvUGSQlYcFDD5U8gxBiKLHvS\/eFE489nrSdJ5DIM6jP\/firwqkHH968tc1BUhxk3J301Ss+NHWqk2Ov+hG3OIiiYVJZIdK8MCl+NoTdSdNCJG9my\/Ot4B4CBAgQIECgDwKCpD5U0RwIEKhMQJCUn7btj7VNm1mWR92qfsStSJhUdog0HibtW10NJ06e3MLX5zBJiJT\/O8CdBAgQIECAQD8FBEn9rKtZESBQkkAfgqRIsfPqr4T1E6cSlTp2JHXpsbZpS+WWN6yE+2++9YyPFu1MuuiCc0taeU82k3VnUlUh0niYdGBtLayvr\/c+TBIilb6cNUiAAAECBAj0QECQ1IMimgIBAtUJ9DFIqiNM6tpjbdNWUNyddOeHrtvyeF68bjxMquMRt9jntDBpWh2rDpHGnXYuLfU6TBIiVfe9qmUCBAgQIECg2wKCpG7Xz+gJEKhYQJCUD7iLj7XNmunk7qrRdaNAqY5H3NKESXWGSCODvoZJQqR8v\/fuIkCAAAECBIYhIEgaRp3NkgCBnAJVB0lxWON9VHHYduxj\/NG2+OcqH2+bFry0\/aDtRctj1kHcTYRJ4wenj8Z9\/fV7z5hCXea7lpfDocOHt\/Tf5TOThEiLfht8ToAAAQIECAxdQJA09BVg\/gQIzBXoS5C062N3hEP3fGdzrlUGSV\/c\/eZw8vY7NvuqK9CoeinPOog7zq+uXUlxjpNv4WsyRBqZX3vttWHvykrnwyQhUtW\/RdonQIAAAXBTkxQAACAASURBVAIE+iAgSOpDFc2BAIHKBPoSJO2+4c5wcP2RWoKkrh+0vWgxzdpx1USY9MGPXnzGcJsK7kZh0o4dO8LRo0eTcXVpZ5IQadHK9zkBAgQIECBA4AkBQZKVQIAAgTkCdQdJcShVPN42eVhzlTuS+h4kxRpN7rqKf3fsk5\/espKqeIPbeAfTzkS68ML9lT62uOjLYjxMitd2JVASIi2qrM8JECBAgAABAk8KCJKsBgIECAiSSl0DQwiSItgtb1gJ99986xl244FSVWHSrBBpNJgqg8JFiyWGSftWV8O27duTS9seJgmRFlXU5wQIECBAgACBrQKCJCuCAAECgqRS18DnfvxV4dSDD2+22dSjVqVOakZjsw7hHoVJVQRJ00KkK97+22H9xKkto2w6TDqwthY2NjZaHSYJker4LdEHAQIECBAg0DcBQVLfKmo+BAiUKlDXo23xXJn4E3dvdP3RtvGgo8oQaVqgMq\/4VY2lzjBp2pxH85p8M1+0aDJMiv3vXFpqbZgkRCr1q1JjBAgQIECAwIAEBEkDKrapEiDQXoEYWI0OKRYkza9T1gBpWmtlh0oxTLrzQ9dt2YkV+407k8ralTQvRBrNUZiU7ndciJTOyVUECBAgQIAAgWkCgiTrggABAi0QqDpIilM8+4rDmzOtcqdKVTuSygiQqgyVqgyT0oRIbQ6Tdi0vh2PHj2\/yx513Tb3RTYjUgi88QyBAgAABAgQ6LSBI6nT5DJ4Agb4I1BEkveB9t4UTjz2ekHUtSFoUIs3aYTTrsbPJdVPWDqXY3zd+78bw6N3f3NJFkfazhEix0\/iGvgOHjrXqzKRkXNdeG66+5ppGwyQhUl++Mc2DAAECBAgQaFJAkNSkvr4JECDw9wKTZzFV8Xjb+GNPVQVJk8FNkQBltDhmhUhZ2571lrXxRZi1zVkLePLNdfG6PG1nDZFG44lh0r4v3bcZHFYdHqb9RY5h0t6VleQxzvhT186kaQFS7P8XH7g97dBdR4AAAQIECBAgMPq3y+mm9pYrAQECBAhsCoyCpCrPSZo8P6eKMKmuIClPKDO+3L64+83h5O13zFyBRduPDU\/rI0u7eUMkYdLWsgqRfNESIECAAAECBMoVsCOpXE+tESBAIJdAHTuSdn3sjnDonu9sjq+KIGly10+W4GQaXNEwZVExpu0cGr+n6Pin7YJK02ZZ8447k\/YevHcLQxV1X+Q8+XncmbRvdTVs2749+SjuTIo\/Zf\/ftoRIWSvjegIECBAgQIDAYgFB0mIjVxAgQKBygTqCpMlQoYpAYXIXTprQZB7uZKBStL1ZfX3ux191xhvXRtcW6TOav\/Ta3wr333zrlq7ntVlWiDTqcKhhkhCp8q8tHRAgQIAAAQIDFRAkDbTwpk2AQLsE+hIkTe7wKRLClB2oLKr4ooO5884lBjnnvf\/dZzxKN629qubc5jDpwNpa2NjYSMpT1s4kh2ovWu0+J0CAAAECBAjkFxAk5bdzJwECBEoVGA+Tqjhsu2s7kurajTRZxHmBUp4wKbrHn0VhUlUh0mh+bQ2T4vh2Li2VFiYJkUr9WtIYAQIECBAgQOAMAUGSRUGAAIGWCPQhSCrrjKSyH5HLWuIyw6RRkBTH8MLLLgmP3v3NLcOJ4VTVIVIXwqRdy8vh2PHjmzZ53ugmRMq60l1PgAABAgQIEMguIEjKbuYOAgQIVCJQ9eNtk6+Er+KMpLLe2lbmI3JFijXvMO60u5PGg6RZYdLkGNO2nWduu2+4Mxxcf2TLrVWshVxj27Mn3HXkSHJr1sfchEh5xN1DgAABAgQIEMguIEjKbuYOAgQIVCIwCpJ27NiR\/CO6isfbdl79lbB+4lQy\/irCg7KCpKYea5tW2DJ2J2UJk6oMkUbzm3yDX1XrIc8vSnyj296VlS23LnqbmxApj7R7CBAgQIAAAQL5BARJ+dzcRYAAgdIF6g6SqggPqgiS6ghWFhWzaJg0GSTtuut\/hq9e8aEzuq1zrn0Jk4RIi1avzwnULxDD3\/vuuy\/8+Z\/\/efjc5z4Xbr\/99vB3f\/d34SMf+Uh4xSteccaAHnjggeSzG2+8MZx99tlh\/\/794fzzz69\/4HokQIAAgVQCgqRUTC4iQGCoApOPmy3aGVHEKfYVdyPFn6p2JE2GB2XvSiojSGr6fKRZNYxz+8bv3XjGGUfx+kUB0GSQ9NzX\/j8zl8qitoqsscl7x3eojT4re03kHe+0nUmxrfHfQSFSXl33EahO4Atf+EK46qqrNjs4cuRI+N73vpf8+WUve1n42Mc+Fs4555zNz7\/+9a+Ht73tbZuPs8YPfu7nfi685z3vCc94xjOqG6iWCRAgQCC3gCApN50bCRAYgkDdQVI0rfLRtqrf3FZGkNSW85Fmre9Z5yYtCoBGYdK8EClNKFX2791kmNSWICnOM4ZJ+1ZXw4mTJ7dMO4ZJQqSyV4L2CJQj8Nhjj4Uf\/MEfDE996lOT4Pf6668PV1xxxWbjn\/rUp8IFF1yQ\/Pmuu+4Kb3nLW8JLX\/rS5M8HDx4M3\/\/+98PFF18c3vrWt4anPOUp5QxKKwQIECBQqoAgqVROjREg0DeBOoOkaNf1N7eVESS16XykWet5ctfU6Lp5YdKsIOlll198xmNui0KpMn\/P4rgOHDq2eXZWbLttYdKBtbWwvr6eTPv65+6cOv1ffOD2Mlm0RYBASQJf+9rXwutf\/\/pw\/O\/fyvgbv\/Eb4Q1veEN45JFHwjvf+c7wohe9KLzjHe8IT3\/608O9994b7r\/\/\/iRYetrTnlbSCDRDgAABAmULCJLKFtUeAQK9EhAkZStn2UFSnYFKtpmGMOvcpFljjoHNtN1I8fppbdU598k3+rUtTIrj2bm0FC5\/6GwhUtaF6noCDQs8\/PDDyaNrt9xySzKS173udSGGSR\/\/+MfDHXfcEd797neHZz\/72Q2PUvcECBAgkEVAkJRFy7UECAxOQJCUveTjO4ryhCFF788+4vx3ZAmTJndaxV7HfYRJ8+swzS\/eYSdS\/vXrTgJ1CMRH1a688srwiU98IukuPtb2pje9KXzgAx8I733ve5Nzkxb9fPe73w2f\/\/znwx\/8wR+EW2+9NTlz6fnPf374hV\/4hXDhhRduOXNpUVs+J0CAAIHiAoKk4oZaIECgxwJNBkmR9dsPPVSq7uTOkyoeYSoaBBW9v1SwFI2lCZMWhUijbtoQJu09eO+WWVexRlKwbrlkmt+Fx\/6\/zWuqPAQ\/61hdT4DAmQJ\/+Id\/mDzGFn\/iuUfxDKVLL700vPa1r114DtLf\/u3fhl\/91V8NX\/7yl6fSPvOZz0wO9\/7pn\/5p9AQIECBQk4AgqSZo3RAg0E2BJoKk0WHbVQRJsc3xw5WrCAmKBkFF729ipc0Lk6aFIMc++elkmBddcO4Zw236sPHJA9mTdXj5TzTBmvQ51e9XXxf2rqxsHkwfrxMmNVYiHRNYKHDbbbclO4ceffTR5Nqf\/MmfDB\/+8IfD8573vLn3PvDAA8nB23\/zN38T\/tW\/+lfhx37sx8Ljjz8e\/uIv\/iLceOONm+cuxbfA\/e7v\/m74Z\/\/sny0ciwsIECBAoLiAIKm4oRYIEOixQBNB0jhn2TuSJoOkKkKCIkFQGWcsNbUcZx3APTmeUYg0K0hqeldSHFdbwqR5O7niG91GYVIc89GjR4VJTS1+\/RJYIPCtb30rvPGNbwx\/9Vd\/lVx52WWXJQduT\/537GQz8Y1vcTdTfATuH\/\/jf7zl43vuuSfZ5RRDpfgTH5eLO5fibic\/BAgQIFCtgCCpWl+tEyDQcYE+Bkm7b7gzHFx\/ZLMyZe82+dyPvyqcevDhzfaznJNUVZA063yd0SCzjHHekp7cTTR5bexn9Pa2WUFS\/Hth0vSdSJN1mhYmRT+7kzr+xWv4vROIO5Hio2x\/9Ed\/lMwtHrj9rne9a+5jbSdPngyXXHJJ8sa3l7\/85VNN4hvhYkAVdy6df\/75YW1tLWzbtq13fiZEgACBtgkIktpWEeMhQKBVAnUHSXHy431WsSNpcrdJn4OkRQHStMVWNFSaFSaN2h0PkuaFSdN2OBUdW9ZfrsnQMd5f9nqZNqa0Z0rFe2OYtG91NWzbvj1pKu5MEiZlrbTrCVQvEHcVffSjH006igduf+QjH5kb+sTH4T7zmc+Eyy+\/PDzjGc+YOsDHHnssCaj+y3\/5L+GHf\/iHw3XXXRd+6Id+qPrJ6IEAAQIDFxAkDXwBmD4BAvMFBEnZV0iRXUVF7p0caZ4QadRGkcBm1nlJse3JHUnzgqT42bRQqsjYslczhF0fuyMcuuc7W26tMkzKEiKNBiVMylNZ9xCoT2B9fT1cdNFFm0Hvueeem4Q+55133sxBfP3rXw\/xbW0\/+qM\/Onego4Dqn\/7Tf5oEy\/\/oH\/2j+iamJwIECAxUQJA00MKbNgEC6QQESemcxq8qEgYVuXd8DJOP180KiNI8ipZVYF6AlTVIasMjbnH+4we0xz9XESTNcksbnMV\/QB5YWwsbGxtJyexMyrpyXU+gGoFvf\/vbya6ieJbRgw8+uHngdnwM7ZWvfGXhTkdB0s\/8zM8kb287++yzC7epAQIECBCYLyBIskIIECAwR6DpICn5R\/tDD5Vao6ofbSsSBhW5dxwp666W2O83fu\/G8Ojd3zzDOm2QEW9Mswsq7TlJo4EMIUwqGiKNF23n0pIwqdRvDI0RyC8Qzyv77Gc\/m5xd9Gu\/9mvJo21\/+Zd\/mTS4d+\/e8La3vS1\/4yEkO5Z+\/dd\/PTmQe2VlJbz5zW9eeIB3oQ7dTIAAAQKJgCDJQiBAgMDAgqQ43Re877Zw4rHHk5mXvbukSBhU5N5RGbOGSOPlj\/3f+aHrthwWHj9PEyZN6\/dll18cvnrFh85YYYve3DZ5QxsecavqTW5lhkgjt13Ly+HY8eObjN7o5mueQDMCd911V3jLW96SBDyvfvWrw5VXXhk+8YlPJIOJu5E+8IEPFNpBdOzYseSRuYceeih5rO3FL35xMxPVKwECBAYmIEgaWMFNlwCBbAJ93JEUBcYfVSo7SIrtj4cDaUKYUVWqCJKy9D8aR56DricDkVG\/sx6fG4VJF11w7sJFOS3gyjOvhR0tuKDsMKlI6LdoLrv37Al3HTkiTFoE5XMCFQmMHmnbvn178va1eGB2PDw77kyKPzt27Agf\/\/jHw4te9KLkz3H30uHDh8NLXvKS8OxnPzvVqG6++ebwS7\/0S+Ed73hH8p9PecpTUt3nIgIECBAoJiBIKubnbgIEei7Q1yBp8gDlssOkrgdJcVnf8oaVcP\/Nt25Z4bPCm1khUrx51uHbWYKkWe00ESaVdV5SlSHSqGhxh8LelZXkH6zxx86knn9hm15rBGIo9Pu\/\/\/vJY2379+\/ffJNafBPbhRdeuHlO0qc+9ankDW7xJ372wQ9+MLz\/\/e8Pz3ve8xbOZfTGtocffji8733vC+ecc87Ce1xAgAABAuUICJLKcdQKAQI9FehrkFT1OUnjh11nDTvyhlDTwpasfU8u4zTnE6UJROaFSWl2JI3GlWY8Vf8qlrErKY1ZWfMYD5McwF2WqnYIPCkQzyk6dOhQeOpTnxp27twZnva0p4VbbrklXHrppeFd73pXiIdgj\/679IEHHghvetObwle\/+tWkgdE5Sffff39yxlHcVfTyl788Fe8f\/\/Efh9\/8zd9MHo87\/\/zzU93jIgIECBAoR0CQVI6jVggQ6KlAE0FSpBzvt+zDtmP7VQdJ449zZQ1zitw7uYsoa9\/TlvG88CZLIDLrEbesY5xsJ+v9ZfyqFgmTspiVMdZkvY\/tTBqFSfHv464JPwQIFBOIv1\/vec97kkae85znhJ\/6qZ9KgqRf\/uVfDq95zWu2PG526tSpJFyKu5XizzOf+czk\/KQ\/+7M\/Cz\/7sz97xvWzRhbPXnr729+enL\/0qle9ygHbxUrobgIECGQWECRlJnMDAQJDEhAk5av25BlDWcKOIkFJlnvPvuLwGZOb9YjfrDBp3iNtk43P2pUUr8vi04ZdSXHMk49Hxr+b94hkFYdqZ1md8R+7+1ZXw7bt25Pb7E7KoudaArMFrrrqqnDgwIEtF8RdR\/Hcongu0uTPn\/zJn4Rf+ZVfCd\/\/\/vc3P5p3\/eT9o91L8fDun\/\/5nw8\/8AM\/oDwECBAgULOAIKlmcN0RINAtgTYESck\/0B96qHS48SCl7DOSihyaXWRXUdog6YwQaf\/yk757\/3RqIDLtAO7xoqQJg6adu5Q1SIrXF\/EtayHFXUn7vnTf5tv\/5gVJTYdIoznHMOnA2lrY2NgQJpW1ELQzeIF77rkn2ZH0+c9\/Pnlr2ute97rkDW3TQqSIFR+F+53f+Z2wtrYWnvvc54Y3vvGNc68fB44HeP\/7f\/\/vwz\/5J\/8kvP71r3e49uBXHwACBJoSECQ1Ja9fAgQ6IdBkkDR+QPCQgqQiIUmanVBbQqTxAGnKipx0n\/V4WtowaNrb10bdpgmixoeYNjSr8hctzSNubQmRxh12Li0Jk6pcGNomUIFAPFz7wx\/+cNi2bVsSIsWzmPwQIECAQDMCgqRm3PVKgEBHBPocJL3gfbdt7iYpe0dSLG\/eQ7NrDZLiQDOESbMeT8sSApXRRhx2Eacyf\/1233BnOLj+yJYmR+upjSHSaKDCpDJXgbYIVCsQH4O7+uqrQ3xD2zvf+c5w9tlnT+0w7li6++67w4\/8yI9UOyCtEyBAYOACgqSBLwDTJ0BgvkCTQVIcWdyVFM9yqWJHUlmvcZ8lWOTQ7KpCqMk5LwqRRnMb989yLtI0m7iL54WXXRIevfubZ3ycJZCKN4+\/HS\/+Oev9Zf3+n+EaQrj++r2F51fW+Ga1s2t5ORw7fnzzvKR4nQO4q1bXPoFsAjFEimcwPfjgg8m5S8961rOmNhB\/d2+88cbk81e+8pXZOnE1AQIECGQSECRl4nIxAQIE6hEYBVhVBklVv7ktzWNmszTrCpL+j+v2hBMnT6Yq6ihMKvrWsei+667\/Gb56xYem9pslDGrLwduTa6kLIdIIf\/RGt\/FiCJNS\/Uq4iEDlAvE8sz\/4gz9IHmm75JJLwjnnnDO1z\/g2uD\/+4z8O8SDuePj3rOsqH7AOCBAgMBABQdJACm2aBAh0T2B8N1QVO5KqDpKKPHqV9\/yfRX1O7pzZ\/w\/\/V9i7spJqccQalPGoVnSPP33blTRaT10KkYRJqZa+iwg0IhAD3c9+9rPhsssu2\/KGt3mD+Xf\/7t+Ff\/2v\/3WY3E3cyAR0SoAAgR4LCJJ6XFxTI0Cg2wJDDpKq2s00+cr6GCTFt3itr68vXCzXP3fn3GvS7iQaBUl925VURsi2sAgVXjBtZ1Lszu6kCtE1TWCOwBe\/+MXw9re\/PTzyyNYz2Gbdcu6554brrrsunHfeeVwJECBAoGIBQVLFwJonQIBAXoGqg6Q4rj4euD1rN1MMCsL\/+aqw9+ATO4LiTwyS4k+aXUmTQdLLLr94y+NpWYOk2O+5b70onHrw4TOWSNq24o1teLxtVoh04YX7QxUHuef9nVp0X1wj+1ZXz3jcUZi0SM7nBAgQIECAwJAEBElDqra5EiDQKYHJrflVPN7W1gO3Fz2iNq+Qs+6dFyT9p898Jhw6fHhms5MhUgx68gY4ox1JsbPz3v\/ucPL2O6b2myVMuuUNK+H+m2\/d0k6W+4v8YswLkUbtdi1MmrZLTZhUZJW4lwABAgQIEOiTgCCpT9U0FwIEeiVQR5BU9TlJRR5RG38jWZZQZG6QFEI48N0LwvqJU5trJc2upGlBUmwgz1lO40FS1x9vm3Xw+OQjhF0KkmJdY+goTOrV16nJECBAgAABAiUKCJJKxNQUAQIEyhToQ5BUZGdRFSHUvF1JsXbTHnGbFSLF67OOcTxEivdfdMG5Z4RR42soS4CWdSxF1+q8t9fFee770n3hxGOPb3bTtTApDnzn0tIZ52fZmVR05bifAAECBAgQ6LqAIKnrFTR+AgR6KxCDpB07doSjR48mc6zi0baqdyQVCZKK3DsrVEmCpBgYfWtpy7oZ7Uqadj7OvCAp6xinBUnTHpEbDS5LkJT3Ubusv0BpD9WeXFvJGr78J7J21\/j1u5aXz3jsUZjUeFkMgAABAgQIEGhQQJDUIL6uCRAgsEig6gO3qw6S4vzGg4cswUjWe+c9Mnbsk59+gvorn3viPycO3Y5\/Ne0Rt3kh0qh2k8HKrDlOC5FiG2UFSbGtqnclpQ2RRja7b7gzHFx\/8o1LXQyS4lx279kTDt50UxLsxp8Y7gqTFn17+ZwAAQIECBDoq4Agqa+VNS8CBHohML4rqYodSRGprQdux7EtOoNoMpwZL\/oLL7skPHr3Nzf\/ajNM+vu\/OXDo2JazkkZh0vjZOHUESdPmOT6PLOFb1h1SWX5JsoZIo7bH3wwY\/66rYVLczRYffRzfJShMyrKCXEuAAAECBAj0RUCQ1JdKmgcBAr0UqCNIqvpg5CK7ZKYFI\/PCo\/FFMHmQ9WSQFK+dPMdnchFdf\/3eLX81LdRJsyNp1m6kUePT3ro2+ixLkJQmfMvzizLvPKRF7dWx623RGMr6XJhUlqR2CBAgQIAAgS4LCJK6XD1jJ0Cg9wJDP3B7MkiaFgbNWwSLdiXFe6ftTBq1OR4kzQp0ygiS5j3eFseSJUwqe1dSkRBp5Fj1rrc6vwhGYdJkn3Yn1VkFfREgQIAAAQJNCgiSmtTXNwECBBYI9CFIilPMe05S3M1y7lsvCqcefHhTalqYFN9+Nu0n7Y6maWFSmt1Ik3ObFvos2o0U72ljkJT3UbZpdejTrqQ4vxgmTTuYXZjkK50AAQIECBAYgoAgaQhVNkcCBDor0JcgadFZR5MFSnNw9qzwaLKteX1PhjzjgVKaIGmy7TYGSdPGtOgXoswQadRXX85KGs0nhknj52mN\/l6YtGh1+ZwAAQIECBDouoAgqesVNH4CBCoVmAxy6v5HYh1BUgSs+tGjtI9bTTv\/aPKsoyyPecW5Tdvtk6aNNI+sLbomzW6k0QKe3LWVdxdXbG\/RuOb90lQRIsX++vIGt0m7nUtLYX19fctf1\/09UemXoMYJECBAgAABAhMCgiRLggABAnME2hIkjd4UVdWb26o+cHtRkLToAO3Js47SBEHjZc0TJqUJYxZdU1aQFOeSZc5Zd4BNC7NGf5el33lfJpOPt8Vru\/oGt8l57lpeDocOHxYm+W8TAgQIECBAYBACgqRBlNkkCRDIK9B0kBTHXceb2+o4w2ZWuDEvRBo9vrYoiEpT32lvRpsXkiwKiRYdQp0lRIrj\/9yPv2rzLKg4rkX9z5tzmkfuJu9fNJ80xouuqXrn26L+q\/w8hknHjh9Pujh69Gjyn3YmVSmubQIECBAgQKApAUFSU\/L6JUCgEwKCpPLKlOUNbNPOP8q7y2Z8BlkClnlBTprQJWuQNG1s42PPsjMoyw6sqh5lm7Zy+rwrKc43npt09TXXbAZJwqTyvj+0RIAAAQIECLRHQJDUnloYCQECLRQQJJVXlLRBUpY3sOUZ3azAZjKomRUkpQlesoZIcR7jO5Linyd3JVURJKWZSx7jWff0PUgahUl7V1a2ENiZVOYq0hYBAgQIECDQtIAgqekK6J8AgVYLtClIilDxkZmqzkkaf+yoirNrYogwftbRsU9+ekvt07yFbVrYkmcBfXH3m8PJ2++Yeut4YDMraJm8cTLkyRMkTTtcO++B22l2JNUdIo3MJh9vi39fxXrLsy7KuifuTJoMk2LbAqWyhLVDgAABAgQINCkgSGpSX98ECLReYEhBUpUHbo+Clck3sMUwKU2ANFooZZyVNGpr2plJeRZkGSFS7HdRkBSvybIrqeijeXks0twzhF1J0SGGSftWV8O27ds96pZmYbiGAAECBAgQ6IyAIKkzpTJQAgSaEGhLkDQ+96p2JFVx4PbkzpzJIClLMBIN0uy0ybJOYnt3fui6zUOus9w7LdiZdnB4mqBscofUyKXMA7enHeA9mm\/WOmR1mrz+Be+7LZx47PHNv+7bjqTRxGKYdGBtLWxsbCR\/5RDuoivH\/QQIECBAgEAbBARJbaiCMRAg0FoBQVL+0sx6G9u5b71oy9vJsvYwK3TJ2s749dMCqnntTQte8oZIsZ9ZB4mXGSTNmk\/dIVIcR5W734qsg6ru3bm0JEyqCle7BAgQIECAQO0CgqTayXVIgECXBNoQJEWv8XG0fUfSrAApziPuzikaBJW9K2lyPaY9P2n8viIhUmxn1tlPRYKkNI\/uNREixfkO5fG28TUiTOrSN7+xEiBAgAABAvMEBEnWBwECBOYIDClIigxFDtxeFCCNmMs456iMNspa+EVDpDiOtG+IyxL8zHo7XewvSztlOU22M3nodl8fbxuf967l5XDs+HFnJlW1qLRLgAABAgQI1CIgSKqFWScECHRVYGhBUt5HjtKGSKN1ML4DJ2+oMetxsDrXWhkh0rw30RXZkdTUW9nS+ldxJlfavpu8btob3bzNrcmK6JsAAQIECBDIKiBIyirmegIEBiXQxiApFqAtj7dlDZBGi6fo422xnaofcVu00MsIkWIfWd6uljZ0a3uIFOc91CApmfu114a9KytblpgwadFvnM8JECBAgACBtggIktpSCeMgQKCVAoKkn5hZl1khUpq3lJX1aFpZ7WRdfGWFSJO7keI4xsOirDuSZgVIk+1mnW9V1w\/x8baR5bQwKX4mUKpqtWmXAAECBAgQKEtAkFSWpHYIEOilwNCCpFjEs684vFnLaefW5N2FNLlAyno0rax20i7gskKk2N+ioGjR5+NjnhcitTVIGvKupFiTGCbtW10NJ06e3LL8hElpfxtdR4AAAQIECDQhIEhqQl2fBAh0RmCIQdKsA7fLCpBGxS9rN1Gdj7jVGSKlCZpGlotCJEFSe79yhEntrY2RESBAgAABAtMFBElWBgECBOYIDDFImnbgdtkh0oi8rN1EdYRJZYZIaUOiRTuS0gRIbQ2R4riGziXyWQAAIABJREFUviNp9HsQw6QDa2thfX19y7eRnUn+64kAAQIECBBoo4AgqY1VMSYCBFoj0JYgKYKMj6Wqw7an\/eN+\/788d2o90pyFtKiQZe1Kiv1MHuBdZoDSRIi0KGyad6D2ogBqUV3q+nwySIr9Tnucsq7xNN3PzqUlYVLTRdA\/AQIECBAgsFBAkLSQyAUECAxZQJAUwmSQVEaANL6mytqVFNucbKtomFTkQPFZvzfTAqBZb2ObFggtOlD7ljeshPtvvnVL92nf9lb37\/ruG+4MB9cf2dLtkIOkCLFreTkcOvzEOWU7duwIR48edQB33QtTfwQIECBAgMBcAUGSBUKAAIE5Am0NkuKQq9qVFMOTfV+6L5x47PFEZjxIKjtEiu2XuSupzDCp6RApziXPo2tlBnNVfzlMe4yy6j670P7uPXvCXUeObA5VmNSFqhkjAQIECBAYjoAgaTi1NlMCBHIIDC1IGoUn\/+lrD4ZD93xnU6zqXSJlhx9FdiZVdR5Ulp1Ik0t10S6k8es\/9+OvCqcefHjzr9q6GykOcPxg9yQcvfwncvyW9vOWeG7S3pWVZFdS\/BEm9bPOZkWAAAECBLooIEjqYtWMmQCB2gSGEiRNC0\/2Hry3tiCp7F1JceB5zkyqYhdSHEveEClLgDQqVlfORxIkLf4aEyYtNnIFAQIECBAgUL+AIKl+cz0SINAhgbYFSeO7E8p4tG3e7ps6g6S4JKrYSTNtZ1Lsa3KXTp27kKb1v+hXIoZDccyj\/5x1fZdCpDiHF7zvts1HKOOf7Ug6s7LTwqR4lTe6Lfqt8TkBAgQIECBQlYAgqSpZ7RIg0AuBtgVJEXV0AG\/RIGleeBL7GT8nqY5\/5FexKymOe1aYNAp0qtiFlGcnURm\/MF0Lks6+4olDpUc\/gqTpqyCGSftWV8O27duTC+JjbsKkMn5jtEGAAAECBAjkERAk5VFzDwECBBoSiMFW0SAp7e6bJg5CHt+VVObZPpMh1WT5jn3y01v+Ku+h4osOxy5zTpNzyPv4XENLOelWkJReP4ZJB9bWwsbGhjApPZsrCRAgQIAAgQoEBEkVoGqSAAECVQkUDZKy7L6J19b9eFtVu5JG9Zh2btJkrbKGPYvCo9h+1jazrp8uhkiCpKxVfuL6nUtLwqR8dO4iQIAAAQIEShIQJJUEqRkCBAjUIZA3SEq7C2l8Dn0MkkYO573\/3eHk7XdUXrKqA6TRBLoYJE2urzgXj7alW5LCpHROriJAgAABAgSqERAkVeOqVQIECFQiMAqSYuPxnJQ05yRl2YU0OegmXs9+yxtWwv0337o5lDLCmFkGL7zskvDo3d8svVZljDntoLoYIsW5CZLSVnj6dbuWl8Ox48eTD52ZVMzS3QQIECBAgEA2AUFSNi9XEyBAoFGB0eHfac5JyrMLaXJyTe1K+sbv3bgl4MkbzCw6UHx0FtK8A7mzFDzvOLP0MX5tV0OkOIfdN9wZDq4\/smXqdiRlWwnx3KSrr7lm86YYKHmbWzZDVxMgQIAAAQLZBQRJ2c3cQYAAgcYEJt8iN2tHUpFdSOOTayJIiv1POxw7S0iTNkCaV8gYLn3vwYfDqQcfnnpZlvFUtWC6HCRNHuYejQRJ2VdKDJP2rqwkh\/DHH2FSdkN3ECBAgAABAtkEBEnZvFxNgACBRgUWBUll7EKanGATj7flDZPKCJAaLXCGzrscIsVpvuB9t4UTjz2+OWMhUobiT1w6LUyKl9idlN\/UneUJjP57y3osz1RLBAgQaFpAkNR0BfRPgACBjALjYdL4jqSydiFNDqepXUlxHJNvWZu1C2hIAVJ06XqIFOdw9hWHtyw1QVLGL4IpYdK+1dWwbfv25BPnJhXzdPeTApP\/B4y8NoKkvHLuI0CAQPsEBEntq4kRESBAYK7AZJBUxS6k8QE0GSTFcUyeXzQeJg0tQBIi+XKYJxB3Jh1YWwsbGxvCJEsll0BZodGszoVJucriJgIECLROQJDUupIYEAECBOYLjP8P\/f0Hvz714tEh0mVZjj+GVPfOkWnnJR375KfnTq3s+ZflWLSdaTuRYpttOK8py9zsRsqile3aLGHS+HeJf+Bnc67j6qKhTp6aFu1zkUueMS1q0+cECBAgUL+AIKl+cz0SIECgkMDkm9smw6QqQpTJg5HbGiZVMfdCxSrxZiFSiZgDaGrn0lKyM2n0iFuccvxHfJqgwD\/261kgaWpRZCR56lhkTHn6KzI\/9xIgQIBAcwKCpObs9UyAAIFcAvF\/6I+\/oWkUJFUZojT1eNv4o2svvfa3wv0337ppNr4rqcq55ypSyTf1NUSKTHWHkiWXptXN7VpeDocObz2LKu2Ahx4KZA1U8nhl7SNt7UbX5RlT1j5cT4AAAQLDFBAkDbPuZk2AQAcFznrnF54Y9f7lM4KkqoOUuoOkWWcfvfCyS8Kjd39zs3pde6Qrz7LrS4gU5+6RtjwroNg9u\/fsCQdvuilXI30KItoY2pQ1pj7VKddCdRMBAgQI1C4gSKqdXIcECHRJYPJ\/6Df1P9g3Q6QpQdL4m9uqtN159VfC+olTm12UvZNk0cHZseNdd\/3PcOeHrgunHny492FSnwKkaSFS\/Luy11CV67\/Lbcdzk\/aurOSaQhPfeYsCljxjWtRmLpyxm4qMKc+9RcfrfgIECBAgUERAkFREz70ECPReoOkgaUuANNLev5z8\/+LjbfH8k7qCpCp2JaUJj+Jcx3dcTTt8u287k4RIvf9qqXWCMUiKP3nCpDwhR5tDm7Lh8\/iUPQbtESBAgACBugUESXWL648AgU4JNBUkTQ2Q\/l7u9FWv2HJgbteCpLTh0WSANL5wpoVJ8fOuB0qzAqSuz23ykbY4H7uR6vsqHAVJecKkPEFJl4KkPPOrr3J6IkCAAAEC7RQQJLWzLkZFgEBLBJoIkmaFSDFAGv2Mj6uuICn2Pf54W5YgoIzwaHJJ9ClMmhcgCZFa8mXQgWGMB0bzhrtvdTWcOHky9Yyyhi1tDJJST9aFBAgQIECAwEIBQdJCIhcQIDBkgTqDpEW7kMbrMHpz2+jV3nWFSbtvuDMcXH9kcyjzwqQqwqNpYdI3fu\/GLQdwj67pwu6kRQFSH0OkOKcsIeSQv3\/SzD1teDRqq47H29IGSVkDqjQeriFAgAABAgSqFxAkVW+sBwIEOixQV5CUZhfStCAp\/l2bzknKEh7FsZf1trlb3rAS7r\/51jNWWpvDpEUhUpvHnvZXekiPtGUNdKYZXnTRRWlpN69L22+eAGnUicAnc1ncQIAAAQIEei0gSOp1eU2OAIGiAlUHSVl2IU0GSeN\/rmtHUuxzPByIO0uaCo8maxsfdZt8o1u8pk2BzKLwqG3jLfL708cQKW1ok9etiiCpSIAU5yFEyltN9xEgQIAAgf4KCJL6W1szI0CgBIEqg6Ssu5DaEiSNn5MUx7T\/X567ULqsnUcLOwohHFx+TasedUsTHvUpQBrVaDJIqvJxtjwBTxWhTZr1Oe+aPGOa1d7Z27YlH43e7ph3bIKkvHLuI0CAAAEC\/RUQJPW3tmZGgEAJAlUESXl3ITUZJE3uOtp78N7N4cwKkuoMjyZLPesg7tF1Ve9SShse9TFAinMqYzdSnnAoy698ntCmjDHl6TfLvMoKkEZ9CpKy6LuWAAECBAgMQ0CQNIw6myUBAjkFyg6SiuxCmpxClW9uW\/S42rQgqcngaFp5Y5g06yDuyevLCJayhEd9CpAmw5W931o6oxx5diOVEdrM+7XPE+jEMeW5L+fXT6bbyg6QYudCpEwlcDEBAgQIEBiMgCBpMKU2UQIE8giUFSSVsQupyiBpUXA02feBQ8fC+olTm3+dJyjIU48898x61C1PW2XcU0ZoVcY4xtsoK7SZFiLt\/4f\/K1f4UtaYxufZ1hCoSD0XBUjjYVDat6kJkYpUxL0ECBAgQKD\/AoKk\/tfYDAkQKCBQRpBU5i6k8akU2ZGUNTga9TvadRTvH9+V1OYgadwshkrfe\/DhcOrBhwusiuy3VhUeTQtb8oQlZYQ2s0KkqJVnTPG+8XHlbSN7tbpxR5YAadb3xqyZ2onUjTVglAQIECBAoCkBQVJT8volQKATAkWCpCp2Ic37B+G8N7flDY6SEOCCMw\/Tju3t+9J94cRjj28OqSth0mjAX9z95iRQevTub1ayFheFR2WEN9MGnidwKTKWaQFSHFfcibQZQOZ4rX0lRelBo3kDpGlTH\/9+Ex71YHGYAgECBAgQqElAkFQTtG4IEOimQN4gqapdSGmCpCKh0azgaFr1uroradZKLLpbaVFwNN5vkeBm0W9SFUHSrDanHaodx9e1UHGRaRs+HwVIs8YiCGpDlYyBAAECBAgMQ0CQNIw6myUBAjkFsgZJVe9Cmhck7T\/49VyzzHtIdpd3JZUV5OQJbWKRyup\/WsHzjinr4hEiZRXLd70AKZ+buwgQIECAAIHqBARJ1dlqmQCBHghkCZJq2YX0zi88qbp\/eYtwmiApb2g0q5Rt35VUZWATTfKGNlnGlbePKn\/9hEhV6j7RtgCpemM9ECBAgAABAvkEBEn53NxFgECPBdK82WjLm5DGw50Jl9NXvSK31LzdTZuNjoVJ04KksoOjyclMBknx80WPNWUJUcb7yxOo5O0rbdHyjClt2229TohUbWUESNX6ap0AAQIECBAoLiBIKm6oBQIEeiSQJkTanO7eP50586wBUqrQaFpvE0FS1cHRtCFMhknjhyyXuTTyhDZFg6Q8fZY55za1NStAShMetmkebRzLovAojtkZSG2snDERIECAAIFhCgiShll3syZAYIpAHSFS7sBoynhjWDU+5nlvbRvdnjZYyRKgTNuVVEWYlGVMFni5AkKkcj1HrQmQqnHVKgECBAgQIFCtgCCpWl+tEyDQEYHUIVLKXUhlBkaRcHKH0\/h4d+zYkSgfPXo07H\/\/+0sRzxra7L7hznBw\/ZHNvssIkrKOoZSJa+QMAY+ylb8oBEjlm2qRAAECBAgQqE9AkFSftZ4IEGixQKogaU6IVObU0jwW17Ygqe2HbpdZnyG1JUQqt9oCpHI9tUaAAAECBAg0IyBIasZdrwQItExgYZBUVYi0fznX2SdtC5JiOSd3JS06dLtlS8BwJgSESOUtCQFSeZb\/f3v3HitXVfcNfKFHakO5FqpyeaFcEiRoJQao2qRRgT9ETAhETCWN2ILUFhAKQsP9oi1vS6OUhoKUYrmYGpREMIohRB77RwOSchN4aksxCAX6AtaWpxT60Ddr6xznzDmnZ2bPnpl9+UzSP9rudfv81mnSb9ZeoycCBAgQIECg9wKCpN7XwAwIEOixwIghUpxfu0FS3aXYjctNc4lu45zj6221V9t69UqYU0k93sgZDS9AyggyhDBSgJTmZz+72emJAAECBAgQIJBOQJCUzk0rAgRKJNBUkFRb70iB0k4Co+HI2vnPZG3utSCpmQu3O1m6L9\/xfHji1f\/pH8KppE5qZ9+3ECkbUwFSNo56IUCAAAECBPIpIEjKZ13MigCBLgq0FCTFecUwaZjAqJ1QKO2SW\/3mtrTjNNPOqaRmlPL3jG9ly6YmAqRsHPVCgAABAgQI5FtAkJTv+pgdAQJdEmg5TKqbVy\/Co3qWPAVJcV5OJXVp02YwzM4CpNi9E2XNIQuQmnPyFAECBAgQIFAOAUFSOepoFQQItClQliAp+c\/\/pk1tarTX3Kmk9vy61doppPalawFS7dXSxh57HTK3v0I9ECBAgAABAgQGCwiS7AoCBAiEEARJ2W6DCYueCevefr+\/UydbsvVtpzcBUjt6\/2orQGrfUA8ECBAgQIBAcQUEScWtnZkTIJCxQJowKQ8nDhrn3esTSbEsTiVlvDkz6M5rbO0jCpDaN9QDAQIECBAgUHwBQVLxa2gFBAhkKNBKmJSHECkuPc65\/tWaPARJcV5OJWW4MdvoSoDUBt6\/mwqQ2jfUAwECBAgQIFAeAUFSeWppJQQIZCRQtDCpFiTF5a9fv77ndyTVytB4Kin+uVfcMtqkTXbjNbYmoYZ5zCXa7flpTYAAAQIECJRTQJBUzrpaFQECGQnUh0rxBFJjyJTHU0l5OZEUS+BUUkYbscVunEJqEazhcQFSe35aEyBAgAABAuUWECSVu75WR4BABwTyGCbVv96WpyDJqaQObMCddClAas9bgNSen9YECBAgQIBANQQESdWos1USIJChgCCpNUxhUmteaZ\/2Gltauf98C9twPeTl5GH6FWpJgAABAgQIEMhOQJCUnaWeCBCokEDewqTafGqXbufpVFLcFsKkzv1wOIWU3tYJpPR2WhIgQIAAAQLVFRAkVbf2Vk6AQBsCgqTW8YRJrZsN12Kk8Ci2c7H58N4CpOz2op4IECBAgACB6gkIkqpXcysmQCAjgTyFSXk\/kVQjFya1t\/kESG367bnnTjvwClt7vloTIECAAAEC1RAQJFWjzlZJgEAHBARJ6VCFSa27CZBaN6tv4QRSe35aEyBAgAABAgTqBQRJ9gMBAgTaEMhbmJTXO5IaiU+757\/Dw+s2D\/hjr2INVGomPIotuA39A1wfHsWfi\/hZv379gIedQGrjHz9NCRAgQIAAgcoKCJIqW3oLJ0AgCwFBUnrFL9\/xfHji1f8RJjUQNhMgCY+G33dDBUiNIZIAKf3PrZYECBAgQIAAAUGSPUCAAIE2BfISJsV5FOVEUo18wqJnwrq33xcmhRAESO39IDqB1J6f1gQIECBAgACBZgUESc1KeY4AAQLDCOQtSIrTjK\/wbNm0qRA1q3KY1Ex4FIvoBNLwW1mAVIgfc5MkQIAAAQIESiQgSCpRMS2FAIHeCeQhTKqdSCpakBTnW7UwSYDU\/s+qAKl9Qz0QIECAAAECBNIICJLSqGlDgACBBoG8BEn10yrKiaQ45\/hNbouf+H+lfs2t2fAoejiBNPw\/MSMFSO4\/8s8zAQIECBAgQKCzAoKkzvrqnQCBCgnEMKmX\/4ltDLOKFCTVwqQb\/uuN8PbW\/x1y1xQtXGklOBIejfwPhQBpZCNPECBAgAABAgS6ISBI6oayMQgQINAFgaIHSbUw6aKHX9upVl4DpVaDI+FRcz8UAqTmnDxFgAABAgQIEOiWgCCpW9LGIUCAQBcE6sOkop1IqucZ7lW3RsJeh0ppwiMBUnM\/CAKk5pw8RYAAAQIECBDotoAgqdvixiNAgEAHBcoSJNWI8hYoCY46uHn\/3XUtQBo\/fnz\/YPFbCGufXr4+2vnVG4EAAQIECBAgkH8BQVL+a2SGBAgQaFqgDK+3DbfYL9\/xfFj39rZh71CqtcvqlFLa0Kh+\/lnNpekNUOAHG08g1YdHcVkCpAIX19QJECBAgACBUgkIkkpVToshQKDqAmUOkmq1bfaUUq\/2gvCoNfn6AGmolgKk1jw9TYAAAQIECBDotIAgqdPC+idAgEAXBaoQJNVzTlj0TFj39vtdFB44lNAoPb0AKb2dlgQIECBAgACBXgoIknqpb2wCBAhkLFC1IKnGd9o9\/x3WvrOto6GS0CibzSpAysZRLwQIECBAgACBXgkIknolb1wCBAh0QKCqQdJwlP\/n\/64e8Fdvb\/3fEdUFRiMSpXpAgJSKTSMCBAgQIECAQO4EBEm5K4kJESBAoD2Bsn1zW3saWvdaQIDU6woYnwABAgQIECCQrYAgKVtPvREgQKDnAoKknpfABEIItQBp\/PjxiYdvYbMtCBAgQIAAAQLlEBAklaOOVkGAQAEEYsDTjW+gEiQVYDOUeIqNAVJtqbUgqRs\/AyXmtTQCBAgQIECAQM8FBEk9L4EJECBQdoHGe4s6\/R\/p2njxJEj8z\/uWTZvKTmx9ORAYKkCqP4XU6X2fAwJTIECAAAECBAhUQkCQVIkyWyQBAr0UaAyS4lw6+Z9qF273strVG7s+QGp8fa3Te7162lZMgAABAgQIEOi9gCCp9zUwAwIEKiDQzVNJcaz6e2mcSKrABuvBEgVIPUA3JAECBAgQIEAgBwKCpBwUwRQIEKiGQLfCJEFSNfZTr1YpQOqVvHEJECBAgAABAvkQECTlow5mQYBABQS6FSRFylqY5I6kCmysLiyxFh7FoWp3bzUO28nXNbuwREMQIECAAAECBAg0KSBIahLKYwQIEMhCoFthkiApi2rpQ4BkDxAgQIAAAQIECDQKCJLsCQIECHRRQJDURWxDpRZoDJBiR40XaTuBlJpXQwIECBAgQIBAoQUESYUun8kTIFBEgW6ESfUnkqKRC7eLuFO6P2cBUvfNjUiAAAECBAgQKJqAIKloFTNfAgQKL9DNICliuSep8Fum4wuoD5DiYI33IDl91PESGIAAAQIECBAgUBgBQVJhSmWiBAiUSaDTYVKt\/1og4ERSmXZPdmtpDJAaexYgZWetJwIECBAgQIBAWQQESWWppHUQKKnA+++\/H15++eWwatWq8PTTT4eXXnop\/P3vfw9vvfVW+NjHPhYOP\/zwcOSRR4YvfelL4Ytf\/GL4xCc+kXxjWd4\/gqS8V6jc8xMglbu+VkeAAAECBAgQ6KSAIKmTuvomQCCVQDwFsWbNmrBs2bLwu9\/9LmzevLm\/n7Fjx4YDDzwwxIBp7dq14YMPPhgwxsSJE8Oll14aPvvZz+Y+UOpkmNTYtxNJqbZi6RoJkEpXUgsiQIAAAQIECHRdQJDUdXIDEiAwnEAMkJ555plw4403JieQap8JEyaEGTNmhBgS7bHHHv1\/\/uGHH4a\/\/e1vYcWKFUnoVAuV+vr6wplnnhl+8IMfDHg+b\/KdDJLiWuv7FyTlrfrdnY8AqbveRiNAgAABAgQIlFlAkFTm6lobgQIJbNy4Mdxyyy3hF7\/4Rdi+fXsy84MPPjhcddVVYfLkyeEjH\/nITlezevXqcPHFFw\/4ivJJkyaFBQsWhHHjxuVWopNhkiApt2Xv2sQESF2jNhABAgQIECBAoDICgqTKlNpCCeRTIJ5CeuKJJ5LX0eLpotrna1\/7WrjiiivCJz\/5yaYnvnLlyjBz5swBr8Kdeuqp4dprrw1jxoxpup9uPtjpIKn+27ecSupmZXs7Vi1Aavz2tdqsXKLd2\/oYnQABAgQIECBQZAFBUpGrZ+4ECi4QTx7F19Kuu+66AXcdTZ8+PVx44YVh9OjRLa0w\/ud4yZIlYf78+QPaXX311WHq1Km5vTOpU2FS7DcGCfGzfv36IEhqaTsV8uH6AKlW9\/qFCJAKWVaTJkCAAAECBAjkSkCQlKtymAyB6gjEEGnx4sXJr9qrbHH1aUOkmlz8hrdp06YNeMUthilLly4NhxxySC6BBUm5LEuhJiVAKlS5TJYAAQIECBAgUGgBQVKhy2fyBIopEIOjRYsWJb\/qP1m8hrZ169Zw2WWXhQcffHBA3\/GE06xZsyp1KsmJpGL+fLQyawFSK1qeJUCAAAECBAgQyEJAkJSFoj4IEGhaIIZI9913X7jhhhsGnEQ67rjjws0335zJxdjx0u6FCxcOmNPnPve5cPvtt4d999236bl288FOnEoSJHWzgt0dyyXa3fU2GgECBAgQIECAwH8EBEl2AwECXRX49a9\/nZwYqn+dbffdd09ecYvfspbF54EHHgizZ88e0NVuu+0Wli9fHo455pgshuhIH50Mk9yR1JGSdb1TAVLXyQ1IgAABAgQIECDQICBIsiUIEOiawJo1a8KMGTMG3F8UBz\/vvPOSX319fZnMZdWqVWHKlCmD+vrxj38cvvWtb2UyRic6ESR1QrUcfQqQylFHqyBAgAABAgQIlEFAkFSGKloDgQIIvPPOO+GCCy4IK1euHDDbT3\/608krZwcccEBmqxguSIoXec+ZMye39yRFgBgmZfnNWvXhlG9ty2yLda0jAVLXqA1EgAABAgQIECDQpIAgqUkojxEgkF4gBiNLliwJ8+fPH9TJ3LlzwxlnnJG+8yFaPvroo8m3vzV+TjnllDBv3rwwevToTMfLc2eNp5yESXmu1n\/mJkAqRp3MkgABAgQIECBQRQFBUhWrbs0Euizw3HPPhWnTpoWNGzcOGDneifTTn\/407L333pnOaMWKFcnJI0HSv0441X8ESZlutcw7EyBlTqpDAgQIECBAgACBjAUESRmD6o4AgYECW7duTS7XfvDBBwfRxBApnhLK+jPUt7bFMT7\/+c8nr9FlHVxlPf8s+xMkZanZub4ESJ2z1TMBAgQIECBAgEC2AoKkbD31RoBAg0C8EylesP3uu+8O+JuJEycm39SWdagTvw3u+uuvD3ffffegWpxwwglh4cKFYcyYMZWpkyApv6UeKTyKM8\/yvqz8SpgZAQIECBAgQIBAkQQESUWqlrkSKJjAzk4jXXPNNWHq1KmZr2jLli3hoosuCo888sigvqt4R1JEcOF25tusrQ4FSG3xaUyAAAECBAgQINBjAUFSjwtgeAJlFnjyySfDd77znUGnkfbbb79w1113hfiNbVl\/XnnllXDWWWeFl156aVDX3\/zmN8N1110Xdt1116yHzXV\/gqR8lEeAlI86mAUBAgQIECBAgEB7AoKk9vy0JkBgGIH4itmNN94Yli5dOuiJ008\/PXn9bNSoUZn7rVq1KjnpFMdv\/MycOTPMnj078zHz3qEgqbcVEiD11t\/oBAgQIECAAAEC2QoIkrL11BsBAv8WeP3118P3vve98Oyzzw4ymTt3bjjjjDM6YnXnnXeGG264Yci+Fy1aFE4++eSOjJvnTt2T1JvqjBQguf+oN3UxKgECBAgQIECAQHsCgqT2\/LQmQGAYgcceeyycffbZg04G7b\/\/\/mHZsmXhiCOOyNxuZ3cyjR07Nhn36KOPznzcvHcoSOpuhQRI3fU2GgECBAgQIECAQHcFBEnd9TYagcoI3HTTTcm3sjV+Jk2aFG655Zawxx57ZG7x8ssvh2nTpoX169d3ddzMF5Jxh4KkjEGH6U6A1B1noxAgQIAAAQIECPRWQJDUW3+jEyilwM6+OW369Olhzpw5A75JLCs0uLsAAAAbr0lEQVSEBx98MFxwwQVDdhe\/yW3WrFlZDdXTfmIw1MprUYKkzpZLgNRZX70TIECAAAECBAjkS0CQlK96mA2BUghs2LAhORn04osvDlpPPKl06qmnZr7Obdu2hSuvvDLcf\/\/9g\/ouy2tt9YGQICnzLdRyh7UAafz48UnbxpNwrdSo5cE1IECAAAECBAgQINAjAUFSj+ANS6DMAs8991w466yzwltvvTVgmX19feHee+8Nxx57bObL\/+tf\/5qM+dprrw3q+6STTgoLFiwIY8aMyXzcbnXYeKoojttKUOGb27KrVH2ANNRrlK3UJbtZ6YkAAQIECBAgQIBAdwQESd1xNgqBSgmsWrUqTJkyZdCaO3nR9vLly8M111wzaMwYXv3sZz8LkydPLnwNGsOkVgILQVL75fcKW\/uGeiBAgAABAgQIECi+gCCp+DW0AgK5ExguSJowYUJYunRp2GeffTKd8zvvvBNmzpwZ4riNnxNOOCHE1+l23333TMfsVWdpwyT3JKWvmAApvZ2WBAgQIECAAAEC5RMQJJWvplZEoOcC3Q6ShrtkO4ZHt912W5g4cWLPTbKaQDtBUv1dPls2bcpqSqXtR4BU2tJaGAECBAgQIECAQBsCgqQ28DQlQGBogeGCpJNPPjnMmzcv7LbbbpnRxdNI8ZvaVq5cOajPc889N8Rva4uvt5XpkyZMim0ESc3tAgFSc06eIkCAAAECBAgQqKaAIKmadbdqAh0VGO6y7fia2cKFCzO99HrFihVhzpw5g9Zz3HHHhZtvvjmMGzeuo2vtReeCpM6oC5A646pXAgQIECBAgACBcgkIkspVT6shkAuBV155JfkGtZdeemnAfLK+I+nVV18N55xzTnjhhRcGjBNfaVu8eHGYNGlSLjw6MYlWw6Ta8\/FUUvymMa+2\/acqAqRO7FB9EiBAgAABAgQIlFVAkFTWyloXgR4KbNmyJXml7JFHHhkwiyOPPDK5bPtTn\/pU27Pbvn17crppyZIlA\/qKr7H96Ec\/CqeffnpoDFvaHjRHHbQaJMWp115vEyT9q5ACpBxtaFMhQIAAAQIECBAojIAgqTClMlECxRK45ZZbkqCn\/jN27NiwbNmycPTRR7e9mMceeyzMmjUrvPvuu\/19xRDpkksuSU5Dle1epKHAWg2T6p+v8okkAVLbP346IECAAAECBAgQqLCAIKnCxbd0Ap0UWL16dZg6deqAoCeOt2jRohAv3W7ns2bNmjBjxozkFa3aJwZHM2fOTH5VIUSK6xYktbaLBEiteXmaAAECBAgQIECAwFACgiT7ggCBjghs27YtXHnlleH+++8f0H985ez6668Po0aNSjXum2++Gc4\/\/\/zw+OOPDwiRrrjiijBlypTKhEi1xbcSJjU+W5VTSQKkVD9qGhEgQIAAAQIECBAYUkCQZGMQINAxgfjtbdOmTQsbN27sHyNehH3XXXeFY445puVx4+Xas2fPHhAixdfl5s6dG7761a+W+k6k4bAEScNvIwFSyz9iGhAgQIAAAQIECBAYUUCQNCKRBwgQSCuwY8eO5ETS5ZdfHuLl2LVPfLUthj9jxoxpuuunn346uf9o7dq1\/W0mTpyYXKwdv4msyp9mw6TaZdvRqqwXbo8UHsW1x33pQ4AAAQIECBAgQIBAOgFBUjo3rQgQaFIgBkiLFy9OftWHSaeddlry6tsee+yx057iZdrxgu54t9IHH3yQPBtPIf3whz8M3\/jGN1K\/Itfk9AvxmCBp4Dew1YLF+ju0BEiF2MomSYAAAQIECBAgUAABQVIBimSKBIou8OGHH4bf\/OY34eqrrw6bN2\/uX87BBx8cLr300jBp0qQBp5Pi86+\/\/np46KGHws9\/\/vOwYcOGpE18Le7ss89O7kLaZ599is6S6fxbDZPKciKp\/gSSACnTLaUzAgQIECBAgAABAkMKCJJsDAIEuiYQw6Hbb7893Hvvvf2ni2qDH3DAAWHfffcN7733Xli3bt2A00vHHntsOOOMM8KJJ56YhEk+gwWqFiQJkPwUECBAgAABAgQIEOiNgCCpN+5GJVBpgS1btoTVq1eHRx55JDz11FPh5ZdfTk4q9fX1hcMOOyzstdde4fjjjw\/xDqSjjjpqxNffKo1Zt\/hmwqTaPUlFPZE0UoDk\/iM\/DQQIECBAgAABAgQ6KyBI6qyv3gkQINA1gVaDpDixLZs2dW1+7QwkQGpHT1sCBAgQIECAAAEC2QkIkrKz1BMBAgR6LjBSmFT7+3ifUBFOJQmQer6lTIAAAQIECBAgQIDAAAFBkg1BgACBEgnUB0lDveZVlCBJgFSiTWkpBAgQIECAAAECpRIQJJWqnBZDgACBnQs0nljK26tttQCp9g1scTXx5FTt4w4kO5wAAQIECBAgQIBAbwUESb31NzoBAgS6KpDXE0kCpK5uA4MRIECAAAECBAgQSC0gSEpNpyEBAgSKKZCnb24TIBVzD5k1AQIECBAgQIBAdQUESdWtvZUTIFBRgTwESe5Aqujms2wCBAgQIECAAIHCCwiSCl9CCyBAgEBrAr0MkuoDpKFm7Q6k1mrpaQIECBAgQIAAAQLdFhAkdVvceAQIEOixQC1IitOIF1l348JtAVKPi254AgQIECBAgAABAhkJCJIygtQNAQIEiiLQzW9uEyAVZVeYJwECBAgQIECAAIHmBARJzTl5igABAqUR6EaQJEAqzXaxEAIECBAgQIAAAQIDBARJNgQBAgQqJtDJIEmAVLHNZLkECBAgQIAAAQKVExAkVa7kFkyAAIEQ6sOkLO5IqgVI48ePT3jj3Uv1H5do23UECBAgQIAAAQIEyiEgSCpHHa2CAAECLQlkFSQJkFpi9zABAgQIECBAgACBwgsIkgpfQgsgQIBAawKNr7bF1q2eShIgtWbuaQIECBAgQIAAAQJlERAklaWS1kGAAIERBIYKkIZqsrNQSYBkmxEgQIAAAQIECBCotoAgqdr1t3oCBCok0GyQFEnqw6T6C7TdgVShDWOpBAgQIECAAAECBIYQECTZFgQIEKiAQCsh0lActQAp\/l39Rdou0a7A5rFEAgQIECBAgAABAnUCgiTbgQABAiUXaCdEEiCVfHNYHgECBAgQIECAAIEWBQRJLYJ5nAABAkUTaDdIcgKpaBU3XwIECBAgQIAAAQKdExAkdc5WzwQIEMiFQDtBUlyA19dyUUaTIECAAAECBAgQIJALAUFSLspgEgQIEOiMQLshkiCpM3XRKwECBAgQIECAAIGiCgiSilo58yZAgEATAoKkJpA8QoAAAQIECBAgQIBA0wKCpKapPEiAAIHiCQiSilczMyZAgAABAgQIECCQZwFBUp6rY24ECBDIQKCdMMn9SBkUQBcECBAgQIAAAQIESiQgSCpRMS2FAAECQwkIkuwLAgQIECBAgAABAgSyEhAkZSWpHwIECORUQJCU08KYFgECBAgQIECAAIECCgiSClg0UyZAgECrAmnCJK+1tarseQIECBAgQIAAAQLlFxAklb\/GVkiAAIFEoJUwSYhk0xAgQIAAAQIECBAgMJSAIMm+IECAQIUEmg2TBEkV2hSWSoAAAQIECBAgQKAFAUFSC1geJUCAQJkE6kOlGBw1hkzCpDJV21oIECBAgAABAgQIZCMgSMrGUS8ECBAovIAgqfAltAACBAgQIECAAAECHRcQJHWc2AAECBAojoAwqTi1MlMCBAgQIECAAAECvRAQJPVC3ZgECBDIqYAgKaeFMS0CBAgQIECAAAECOREQJOWkEKZBgACBvAgIk\/JSCfMgQIAAAQIECBAgkD8BQVL+amJGBAgQ6KlA4yXcPZ2MwQkQIECAAAECBAgQyJWAIClX5TAZAgQIECBAgAABAgQIECBAgEB+BQRJ+a2NmREgQIAAAQIECBAgQIAAAQIEciUgSMpVOUyGAAECBAgQIECAAAECBAgQIJBfAUFSfmtjZgQIECBAgAABAgQIECBAgACBXAkIknJVDpMhQIAAAQIECBAgQIAAAQIECORXQJCU39qYGQECBAgQIECAAAECBAgQIEAgVwKCpFyVw2QIECBAgAABAgQIECBAgAABAvkVECTltzZmRoAAAQIECBAgQIAAAQIECBDIlYAgKVflMBkCBAgQIECAAAECBAgQIECAQH4FBEn5rY2ZESBAgAABAgQIECBAgAABAgRyJSBIylU5TIYAAQIECBAgQIAAAQIECBAgkF8BQVJ+a2NmBAgQyK3ALrvsksxtx44duZ2jiREgQIAAAQIECBAgkL2AICl7Uz0SIECgtAK1AKm2QEFSaUttYQQIECBAgAABAgSGFBAk2RgECBAg0JKAMKklLg8TIECAAAECBAgQKJWAIKlU5bQYAgQIdF5AkNR5YyMQIECAAAECBAgQyKuAICmvlTEvAgQI5FhAmJTj4pgaAQIECBAgQIAAgQ4KCJI6iKtrAgQIlFVAkFTWyloXAQIECBAgQIAAgZ0LCJLsEAIECBBIJSBMSsWmEQECBAgQIECAAIFCCwiSCl0+kydAgEDvBARJvbM3MgECBAgQIECAAIFeCQiSeiVvXAIECJRAQJhUgiJaAgECBAgQIECAAIEWBARJLWB5lAABAgQGCgiS7AgCBAgQIECAAAEC1RIQJFWr3lZLgACBzAWESZmT6pAAAQIECBAgQIBAbgUESbktjYkRIECgGAKCpGLUySwJECBAgAABAgQIZCEgSMpCUR8ECBCouIAwqeIbwPIJECBAgAABAgQqIyBIqkypLZQAAQKdExAkdc5WzwQIECBAgAABAgTyJCBIylM1zIUAAQIFFhAmFbh4pk6AAAECBAgQIECgSQFBUpNQHiNAgACBnQsIkuwQAgQIECBAgAABAuUXECSVv8ZWSIAAga4JCJO6Rm0gAgQIECBAgAABAj0RECT1hN2gBAgQKKeAIKmcdbUqAgQIECBAgAABAjUBQZK9QIAAAQKZCgiTMuXUGQECBAgQIECAAIFcCQiSclUOkyFAgEDxBQRJxa+hFRAgQIAAAQIECBAYTkCQZG8QIECAQOYCwqTMSXVIgAABAgQIECBAIBcCgqRclMEkCBAgUC6B+iBpx44d5Vqc1RAgQIAAAQIECBCosIAgqcLFt3QCBAgQIECAAAECBAgQIECAQCsCgqRWtDxLgAABAgQIECBAgAABAgQIEKiwgCCpwsW3dAIECBAgQIAAAQIECBAgQIBAKwKCpFa0PEuAAAECBAgQIECAAAECBAgQqLCAIKnCxbd0AgQIECBAgAABAgQIECBAgEArAoKkVrQ8S4AAAQIECBAgQIAAAQIECBCosIAgqcLFt3QCBAgQIECAAAECBAgQIECAQCsCgqRWtDxLgAABAgQIECBAgAABAgQIEKiwgCCpwsW3dAIECBAgQIAAAQIECBAgQIBAKwKCpFa0PEuAAAECBAgQIECAAAECBAgQqLCAIKnCxbd0AgTKIbBjx47wxhtvhD\/96U\/hoYceCk899VR47733wq233hq+8pWvDFrkxo0bk7\/71a9+FcaMGRMWLlwYjj\/++HJgWAUBAgQIECBAgAABAh0VECR1lFfnBAgQ6KzAo48+Gm666ab+QdauXRs++OCD5Pef+cxnwh133BH222+\/\/r9\/4YUXwqxZs8L69ev7\/+yUU04J8+bNC6NHj+7sZPVOgAABAgQIECBAgEDhBQRJhS+hBRAgUGWBrVu3ho9+9KNh1113DfFk0vLly8O1117bT3LfffeFiRMnJr9fs2ZNmDFjRjjqqKOS3z\/88MNh+\/bt4YILLggzZ84MfX19Vaa0dgIECBAgQIAAAQIEmhAQJDWB5BECBAgUReC5554LZ511VnjrrbeSKV9xxRXhu9\/9bti8eXOYPXt2OOSQQ8KFF14YPv7xj4fXXnstvPnmm0mwNGrUqKIs0TwJECBAgAABAgQIEOihgCCph\/iGJkCAQNYC\/\/znP5NX11auXJl0feaZZyZh0tKlS8Pzzz8frrvuurD33ntnPaz+CBAgQIAAAQIECBCoiIAgqSKFtkwCBKohEF9Vu\/7668Pdd9+dLDi+1jZ9+vSwYMGCcOONNyb3JrXy+fDDD8Pjjz8e7rnnnnDppZeGgw46qJXmniVAgAABAgQIECBAoGQCgqSSFdRyCBAg8MADDySvscVPvPco3qF02WWXhW9\/+9tN34MUA6Snn346zJ8\/P6xatSoceuihYdmyZYIk24sAAQIECBAgQIBAxQUESRXfAJZPgED5BFavXh2mTp0a3n333WRxxx13XLj55pvDuHHjRlxsvLD7L3\/5S3ICaePGjckrcvFb4ARJI9J5gAABAgQIECBAgEAlBARJlSizRRIgUCWBDRs2hGnTpoUXX3wxWfbll1+eXLi9yy67jMiwbdu28MYbb4QDDzwwxFNJtdfkBEkj0nmAAAECBAgQIECAQCUEBEmVKLNFEiBQJYF4Eim+yvbb3\/42WXa8cPuqq65q+rW2eqt4r9Jtt93mRFKVNpC1EiBAgAABAgQIENiJgCDJ9iBAgEAJBWoBUFxavHD71ltvDXvuuWfLKxUktUymAQECBAgQIECAAIFSCwiSSl1eiyNAoIoC69atC+ecc05Yv359svz9998\/uSj7iCOOaJlDkNQymQYECBAgQIAAAQIESi0gSCp1eS2OAIGqCWzZsiVcffXV4fHHHw\/\/+Mc\/+i\/cXrJkSTjppJNa5hAktUymAQECBAgQIECAAIFSCwiSSl1eiyNAoEoC8RvX7r\/\/\/hBDozlz5iR3Gz355JMJwUUXXRRmzZrVMocgqWUyDQgQIECAAAECBAiUWkCQVOryWhwBAlUSWLNmTZgxY0Y499xzw6mnntr\/jWvRIJ5GWrBgQRgzZkxLJIKklrg8TIAAAQIECBAgQKD0AoKk0pfYAgkQqIJA7ZW2vfbaK1x88cVh9OjRYcWKFcnJpPgZP358WLp0aTjkkEOS38fTS3\/+85\/D4YcfHvbee+9hiQRJVdg91kiAAAECBAgQIECgeQFBUvNWniRAgEAuBWIo9Mtf\/jJ5rW3hwoXhoIMOSua5evXqMHXq1P57ku67777kG9xqf\/eTn\/wkzJ8\/P4wbN06QlMvKmhQBAgQIECBAgACB\/AkIkvJXEzMiQIDAsALbtm0LTzzxRNh1113DhAkTwqhRo8LKlSvDZZddFq666qpw4oknhl122SVpv3HjxjB9+vTw7LPPJr+v3ZP05ptvhksuuSScffbZYdKkSTvVdiLJZiRAgAABAgQIECBAoF5AkGQ\/ECBAoEACt99+e5g3b14y47Fjx4bJkycnQdL3v\/\/9MGXKlNDX19e\/mvfffz8Jl+JppfjZfffdk\/uT\/vjHP4aTTz550PNDMQiSCrQ5TJUAAQIECBAgQIBAFwQESV1ANgQBAgSyErjpppvC4sWLB3QXTx1deOGFyb1IjZ\/f\/\/734fzzzw\/bt2\/v\/6udPd\/YXpCUVeX0Q4AAAQIECBAgQKAcAoKkctTRKggQqIjAq6++mpxI+sMf\/hAOO+ywcOaZZybf0DZUiBRJ4qtwd955Z1iyZEnYd999w7Rp03b6vCCpIhvJMgkQIECAAAECBAikFBAkpYTTjAABAlUQcCKpClW2RgIECBAgQIAAAQLNCwiSmrfyJAECBConIEiqXMktmAABAgQIECBAgMBOBQRJNggBAgQIDCsgSLI5CBAgQIAAAQIECBCoFxAk2Q8ECBAgMKTAjh07wty5c8Mdd9wRDj300LBs2bJw0EEH0SJAgAABAgQIECBAoMICgqQKF9\/SCRAgsDOB119\/PZx33nnhySefDH19fWH58uVh4sSJ0AgQIECAAAECBAgQqLCAIKnCxbd0AgQIDCWwatWq5FvhHn744bBhw4b+R8aOHRu+\/vWvhy984Qth8uTJYdSoUQAJECBAgAABAgQIEKiYgCCpYgW3XAIECBAgQIAAAQIECBAgQIBAWgFBUlo57QgQIECAAAECBAgQIECAAAECFRMQJFWs4JZLgAABAgQIECBAgAABAgQIEEgrIEhKK6cdAQIECBAgQIAAAQIECBAgQKBiAoKkihXccgkQIECAAAECBAgQIECAAAECaQUESWnltCNAgAABAgQIECBAgAABAgQIVExAkFSxglsuAQIECBAgQIAAAQIECBAgQCCtgCAprZx2BAgQIECAAAECBAgQIECAAIGKCQiSKlZwyyVAgAABAgQIECBAgAABAgQIpBUQJKWV044AAQIECBAgQIAAAQIECBAgUDEBQVLFCm65BAgQIECAAAECBAgQIECAAIG0AoKktHLaESBAgAABAgQIECBAgAABAgQqJiBIqljBLZcAAQIECBAgQIAAAQIECBAgkFZAkJRWTjsCBAgQIECAAAECBAgQIECAQMUEBEkVK7jlEiBAgAABAgQIECBAgAABAgTSCgiS0sppR4AAAQIECBAgQIAAAQIECBComIAgqWIFt1wCBAgQIECAAAECBAgQIECAQFoBQVJaOe0IECBAgAABAgQIECBAgAABAhUTECRVrOCWS4AAAQIECBAgQIAAAQIECBBIKyBISiunHQECBAgQIECAAAECBAgQIECgYgKCpIoV3HIJECBAgAABAgQIECBAgAABAmkFBElp5bQjQIAAAQIECBAgQIAAAQIECFRMQJBUsYJbLgECBAgQIECAAAECBAgQIEAgrYAgKa2cdgQIECBAgAABAgQIECBAgACBigkIkipWcMslQIAAAQIECBAgQIAAAQIECKQVECSlldOOAAECBAgQIECAAAECBAgQIFAxAUFSxQpuuQQIECBAgAABAgQIECBAgACBtAKCpLRy2hEgQIAAAQIECBAgQIAAAQIEKiYgSKpYwS2XAAECBAgQIECAAAECBAgQIJBWQJCUVk47AgQIECBAgAABAgQIECBAgEDFBARJFSu45RIgQIAAAQIECBAgQIAAAQIE0goIktLKaUeAAAECBAgQIECAAAECBAgQqJiAIKliBbdcAgQIECBAgAABAgQIECBAgEBaAUFSWjntCBAgQIAAAQIECBAgQIAAAQIVExAkVazglkuAAAECBAgQIECAAAECBAgQSCsgSEorpx0BAgQIECBAgAABAgQIECBAoGICgqSKFdxyCRAgQIAAAQIECBAgQIAAAQJpBQRJaeW0I0CAAAECBAgQIECAAAECBAhUTECQVLGCWy4BAgQIECBAgAABAgQIECBAIK2AICmtnHYECBAgQIAAAQIECBAgQIAAgYoJCJIqVnDLJUCAAAECBAgQIECAAAECBAikFRAkpZXTjgABAgQIECBAgAABAgQIECBQMQFBUsUKbrkECBAgQIAAAQIECBAgQIAAgbQCgqS0ctoRIECAAAECBAgQIECAAAECBComIEiqWMEtlwABAgQIECBAgAABAgQIECCQVkCQlFZOOwIECBAgQIAAAQIECBAgQIBAxQQESRUruOUSIECAAAECBAgQIECAAAECBNIK\/H+lu529+v1MVAAAAABJRU5ErkJggg==","height":313,"width":520}}
%---
