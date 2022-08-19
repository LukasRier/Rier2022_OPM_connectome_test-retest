function res=hbf_CheckMesh(input1,varargin)
% HBF_CHECKMESH performs some tests on a boundary mesh
%
% status=HBF_CHECKMESH(mesh,lasttesttodo)
% status=HBF_CHECKMESH(points,elements,lasttesttodo)
%   mesh:   hbf struct for triangle mesh
%   points: mesh vertices, [N x 3]
%   elements:   mesh triangle description, [M x 3]
%   lasttesttodo: last test to do --- between 1 and 8.
%
% Tests 1--4 are basic tests that should be carried out always.
% Tests 5 -- 7 test for a closed mesh. If some of tests 2--4 fail, 5--7
% will probably not work correctly.
%
% test 1: does each vertex belong to some triangle?
% test 2: does each triangle consist of three different vertices?
% test 3: does each triplet of vertices form one and only one triangle?
% test 4: does each vertex have a unique position?
% test 5: does each vertex belong to at least 3 triangles?
% test 6: is the solid angle spanned by the mesh at a point outside the mesh 0?
% test 7: does each triangle-side belong to exactly two triangles?
%
% This tool only spots basic errors; it does not attempt to correct them.
%
% v160229 Matti Stenroos

if isstruct(input1)
    p=input1.p;
    e=input1.e;
else
    p=input1;
    e=varargin{1};
end
if isstruct(input1) && ~isempty(varargin)
    lasttesttodo=varargin{1};
elseif length(varargin)>1
    lasttesttodo=varargin{2};
else
    lasttesttodo=7;
end

res.failedtests=[];
res.success=1;
nop=size(p,1);
fprintf('Testing mesh...');
while 1 %... an elegant eternal loop!
    %test 1: does each vertex belong to some triangle?
    % fprintf('Test 1...');
    orphans=FindOrphans(e,nop);
    if ~isempty(orphans)
        fprintf('\nTest 1 failed.');
        res.success=0;
        res.test1.problemvertices=orphans;
        res.failedtests=[res.failedtests;1];
    end
    if lasttesttodo==1
        break;
    end
    %test 2: does each triangle consist of three different vertices?
    % fprintf('Test 2...');
    problemtri=FindUnTriangles(e);
    if ~isempty(problemtri)
        fprintf('\nTest 2 failed.');
        res.test2.problemtri=problemtri;
        res.success=0;
        res.failedtests=[res.failedtests;2];
    end
    if lasttesttodo==2
        break;
    end
    
    %test 3: does each triplet of vertices form one and only one triangle?
    % fprintf('Test 3...');
    problemtri=FindNonUniqueTriangles(e);
    if ~isempty(problemtri)
        fprintf('\nTest 3 failed.');
        res.test3.problemtri=problemtri;
        res.success=0;
        res.failedtests=[res.failedtests;3];
        
    end
    if lasttesttodo==3
        break;
    end
    
    %test 4: does each vertex have a unique position?
    % fprintf('Test 4...');
    problemvertices=FindNonUniqueVertices(p,e);
    if ~isempty(problemvertices)
        fprintf('\nTest 4 failed.');
        res.test4.problemvertices=problemvertices;
        res.success=0;
        res.failedtests=[res.failedtests;4];
    end
    if lasttesttodo==4
        break;
    end
    
    %test 5: does each vertex belong to at least 3 triangles?
    [ntri,ntri_s,ntri_n]=TrianglesForNodes(e,nop);
    if any(ntri_n<3)
        fprintf('\nTest 5 failed.');
        problemvertices=find(ntri_n<3);
        if ~isempty(problemvertices)
            res.success=0;
            res.test5.problemvertices=orphans;
            res.failedtests=[res.failedtests;5];
        end
    end
    if lasttesttodo==5
        break;
    end
    
    %test 6: does the solid angle spanned by the mesh at a point outside the mesh sum to 0?
    orientation=hbf_CheckTriangleOrientation(p,e);
    %if normals point to wrong direction, just flip them...
    if orientation==0
        fprintf('\nTest 6 failed: normals arbitrary or mesh not simply connected.');
    elseif orientation==2
        fprintf('\nTest 6: wrong normal orientation, run "CorrectNormalOrientation.m"');
    end
    if orientation~=1,
        res.success=0;
        res.failedtests=[res.failedtests;6];
    end
    if lasttesttodo==6
        break;
    end
    
    %test 7: does each triangle-side belong to exactly two triangles?
    [stri1,stri0,stri_more]=hbf_SharedTriangleSides(e);
    test7=any(stri1(:)<1);
    if any(test7),
        res.test7.problemtri=unique(union(stri0,stri_more));
        res.failedtests=[res.failedtests;7];
        res.success=0;
    end
    
    break
end
if res.success==1;
    fprintf('OK.\n');
else
    fprintf('Failed tests: %d\n',res.failedtests);
end

function orphans=FindOrphans(e,nop)
trip=unique(e(:));
orphans=setdiff(1:nop,trip);

function problemtri=FindUnTriangles(e)
e1=e(:,1);
e2=e(:,2);
e3=e(:,3);
test2=(e1==e2 | e1==e3 | e2==e3);
if any(test2)
    problemtri=find(test2);
else
    problemtri=[];
end

function problemtri=FindNonUniqueTriangles(e)
noe=size(e,1);
problemtri={};
problemcount=0;
treatedtri=zeros(noe,1);
for E=1:noe,
    if treatedtri(E)
        continue;
    end
    etest1=e(E,1);
    etest2=e(E,2);
    etest3=e(E,3);
    t1=any(e==etest1,2);
    t2=any(e==etest2,2);
    t3=any(e==etest3,2);
    t1(E)=0;t2(E)=0;t3(E)=0;
    test=t1&t2&t3;
    if any(test)
        tris=find(test);
        problemcount=problemcount+1;
        problemtri{problemcount}=[E;tris];
        treatedtri([E;tris])=1;
    end
end

function problemvertices=FindNonUniqueVertices(p,e)
nop=size(p,1);
p1=p(e(:,1),:);
p2=p(e(:,2),:);
sidelensqr=sum((p1-p2).^2,2);
chardistsqr=min(setdiff(sidelensqr,0));%remove zero in case test 2 has failed...
testdist=1e-6*chardistsqr;

treatedvert=zeros(nop,1);
problemvertices={};
problemcount=0;
nop=size(p,1);
for I=1:nop,
    pdif=ones(nop,1)*p(I,:)-p;
    dist=sum(pdif.^2,2);
    dist(I)=1;
    if any(dist<testdist)
        match=find(dist<testdist);
        problemcount=problemcount+1;
        problemvertices{problemcount}=[I,match];
        treatedvert([I;match])=1;
    end
end

function [tr,sides,n]=TrianglesForNodes(triangles,non)
% function [tr,sides,n]=TrianglesForNodes(triangles,non)
% Finds, which triangles belong to the neighborhood of each node
% tr = list of triangles for each node
% sides = index of the node in the triangle (1, 2, or 3)
% n = number of triangles for each node

maxt=20;%maximum number of triangles for one node
tr=zeros(non,maxt);
sides=zeros(non,maxt);
n=zeros(non,1);
not=size(triangles,1);
% count=zeros(non,1);
e1=triangles(:,1);
e2=triangles(:,2);
e3=triangles(:,3);

for I=1:not
    n(e1(I))=n(e1(I))+1;
    tr(e1(I),n(e1(I)))=I;
    sides(e1(I),n(e1(I)))=1;
    
    n(e2(I))=n(e2(I))+1;
    tr(e2(I),n(e2(I)))=I;
    sides(e2(I),n(e2(I)))=2;
    
    n(e3(I))=n(e3(I))+1;
    tr(e3(I),n(e3(I)))=I;
    sides(e3(I),n(e3(I)))=3;
    
end
nmax=max(n);
tr=tr(:,1:nmax);
sides=sides(:,1:nmax);
