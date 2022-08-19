function newmeshes=hbf_SortNestedMeshes(meshes)
%HBF_SORTNESTEDMESHES sorts nested triangle meshes ("small to large") 
%
% meshes_out=HBF_SORTNESTEDMESHES(meshes_in)
%   meshes_in:  the meshes to be sorted, cell array of hbf structs
%   meshes_out: sorted meshes, cell array of hbf structs
%
% v160229 Matti Stenroos

Nmesh=length(meshes);
prange=zeros(Nmesh,3);
for I=1:Nmesh,
    prange(I,:)=max(meshes{I}.p)-min(meshes{I}.p);
end
orders=zeros(Nmesh,3);
for J=1:3
    [foo,orders(:,J)]=sort(prange(:,J));
end
if ~all(orders(:,1)==orders(:,2)) && ~all(orders(:,1)==orders(:,3)),
    fprintf('Meshes not nested. Not changing the order\n');
    newmeshes=meshes;
    return;
end
newmeshes=cell(Nmesh,1);
order=orders(:,1);
for I=1:Nmesh,
    newmeshes{I}=meshes{order(I)};
end
if any (order-(1:Nmesh)')
    fprintf('Sorted meshes.\n');
else
    fprintf('Meshes were already in correct order.\n');   
end
