// Author: OpenAI Codex.
// C++ adaptation of Shutty's Denser Planar Color Codes source release
// (10.5281/zenodo.22820614), including Google LLC's 2023 Chromobius 6.6.6
// construction. Apache-2.0; see docs/superdense_source_NOTICE.md.
#include "circuit_generator.h"
#include "superdense_488_rules.h"

#include <algorithm>
#include <array>
#include <bit>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>

namespace
{
using Point = std::pair<int,int>;
using Qubits = std::vector<uint32_t>;
struct Face {Qubits support; int color;};
struct Schedule
{
    std::vector<std::vector<Point>> layers; // Pairs of qubit indices, not coordinates.
    Qubits prepare_x, measure_x, input_h, output_h;
};
struct Geometry
{
    std::vector<Point> coordinates;
    size_t data;
    std::vector<Face> faces;
    Qubits routing;
    std::vector<Schedule> schedules;
};
int mod4(int x) {return (x%4+4)%4;}

Geometry geometry666(uint32_t distance)
{
    Geometry g;
    const int width = 2*distance-1;
    std::set<Point> data;
    std::vector<Point> centers;
    std::vector<std::vector<Point>> supports;
    const std::vector<Point> offsets{{-1,0},{0,1},{1,1},{2,0},{1,-1},{0,-1}};
    auto inside = [&](int x,int y)
    {
        return x>=0 && y>=0 && x<width && y<width && 2*y<=3*x && 2*y<=3*(width-x);
    };
    for (int x=-1; x<width; x+=2)
        for (int y=x==-1 ? 1:(x/2)%2; y<width; y+=2)
        {
            // Python floor division is needed for x=-1.
            std::vector<Point> support;
            for (auto [dx,dy] : offsets) if (inside(x+dx,y+dy)) support.emplace_back(x+dx,y+dy);
            if (support.size()<4) continue;
            centers.emplace_back(x,y);
            supports.push_back(support);
            data.insert(support.begin(),support.end());
        }
    g.coordinates.assign(data.begin(),data.end());
    g.data = data.size();
    std::map<Point,uint32_t> at;
    for (auto q : g.coordinates) at[q] = at.size();
    Qubits x_ancillas,z_ancillas;
    for (size_t f=0; f<centers.size(); f++)
    {
        const auto [x,y] = centers[f];
        Qubits support;
        for (auto q : supports[f]) support.push_back(at.at(q));
        g.faces.push_back({support,y%3});
        x_ancillas.push_back(g.coordinates.size());
        at[{x,y}] = g.coordinates.size(); g.coordinates.emplace_back(x,y);
        z_ancillas.push_back(g.coordinates.size());
        at[{x+1,y}] = g.coordinates.size(); g.coordinates.emplace_back(x+1,y);
    }
    Schedule s{.layers=std::vector<std::vector<Point>>(8),.prepare_x=x_ancillas,.measure_x=x_ancillas};
    auto contacts = [&](size_t layer,const Qubits& ancillas,Point control,Point target)
    {
        for (auto a : ancillas)
        {
            auto [x,y] = g.coordinates[a];
            Point c{x+control.first,y+control.second},t{x+target.first,y+target.second};
            if (at.contains(c) && at.contains(t)) s.layers[layer].emplace_back(at.at(c),at.at(t));
        }
    };
    contacts(0,x_ancillas,{0,0},{1,0});
    for (const auto& qs : {x_ancillas,z_ancillas})
    {
        contacts(1,qs,{0,1},{0,0});
        contacts(3,qs,{0,-1},{0,0});
        contacts(4,qs,{0,0},{0,1});
        contacts(6,qs,{0,0},{0,-1});
    }
    contacts(2,x_ancillas,{-1,0},{0,0}); contacts(2,z_ancillas,{1,0},{0,0});
    contacts(5,x_ancillas,{0,0},{-1,0}); contacts(5,z_ancillas,{0,0},{1,0});
    contacts(7,x_ancillas,{0,0},{1,0});
    g.schedules.push_back(std::move(s));
    return g;
}

Geometry geometry488(uint32_t distance)
{
    Geometry g;
    const int m=(distance+1)/4, right=-4*m-1;
    std::vector<int> left;
    for (int j=0; j<4*m-1; j++)
    {
        left.push_back(right-std::min(2*j+3,8*m-2*j-3)+1);
        for (int x=left.back(); x<=right; x++) g.coordinates.emplace_back(x,6+2*j);
    }
    g.data = g.coordinates.size();
    for (int j=0; j<4*m-2; j++)
        for (int x=left[j%2 ? j:j+1]; x<=right+(j%2==0); x++) g.coordinates.emplace_back(x,7+2*j);
    std::map<Point,uint32_t> at;
    for (auto q : g.coordinates) at[q] = at.size();
    struct Entry {std::vector<Point> points; int color; std::array<Point,2> pair;};
    std::vector<Entry> entries;
    const std::vector<Point> square{{0,0},{0,2},{1,0},{1,2}},
        octagon{{0,0},{0,2},{1,0},{1,2},{2,0},{2,2},{3,0},{3,2}},
        upper{{0,2},{1,2},{2,0},{2,2}},lower{{0,0},{1,0},{2,0},{2,2}};
    auto add = [&](const std::vector<Point>& shape,int x,int y,int kind)
    {
        std::vector<Point> points;
        for (auto [dx,dy] : shape)
        {
            Point p{x+dx,y+dy};
            if (!at.contains(p) || at.at(p)>=g.data) return;
            points.push_back(p);
        }
        std::sort(points.begin(),points.end());
        std::array<Point,2> pair{{{x,y+1},{x+1,y+1}}};
        if (kind==1 || kind==2 || (kind==4 && x==right-2 && y==8*m))
            pair = {{{x+1,y+1},{x+2,y+1}}};
        if (kind==0 && ((mod4(x)==3 && x+y==right+4) || (mod4(x)==1 && y-x==-3*right-2)))
            std::swap(pair[0],pair[1]);
        entries.push_back({points,kind==0 ? 0 : (mod4(y)==0 ? 1:2),pair});
    };
    for (int h=0; h<m-1; h++) add(upper,-4*m-5-4*h,6+4*h,3);
    for (size_t q=0; q<g.data; q++)
    {
        const auto [x,y]=g.coordinates[q];
        if (mod4(x-y)==3) add(square,x,y,0);
        if (x%2==0 && mod4(x)==mod4(y)) add(octagon,x,y,2);
    }
    for (int h=0; h<2*m-1; h++) add(square,right-1,6+4*h,1);
    for (int h=0; h<m; h++) add(lower,-8*m+1+4*h,4*m+4+4*h,4);
    auto key = [](const Entry& e)
    {
        int y=e.points[0].second;
        for (auto p:e.points) y=std::min(y,p.second);
        return std::tuple{y,e.points[0].first,e.points.size(),e.points};
    };
    std::sort(entries.begin(),entries.end(),[&](const auto& a,const auto& b){return key(a)<key(b);});
    std::set<uint32_t> paired;
    for (const auto& e : entries)
    {
        Qubits support;
        for (auto q:e.points) support.push_back(at.at(q));
        g.faces.push_back({support,e.color});
        for (auto p:e.pair) paired.insert(at.at(p));
    }
    for (size_t q=g.data; q<g.coordinates.size(); q++) if (!paired.contains(q)) g.routing.push_back(q);
    Schedule s{.layers=std::vector<std::vector<Point>>(12)};
    for (size_t a=g.data; a<g.coordinates.size(); a++)
    {
        auto [x,y]=g.coordinates[a];
        std::string signature;
        for (int dx=-5; dx<=5; dx++) for (int dy=-5; dy<=5; dy++) if (std::abs(dx)+std::abs(dy)<=5)
        {
            const auto p=at.find({x+dx,y+dy});
            signature += p==at.end() ? '0' : p->second<g.data ? '1':'2';
        }
        const auto rule=std::find_if(superdense_488_rules.begin(),superdense_488_rules.end(),[&](const auto& r)
        {return r.x==mod4(x-right) && r.y==mod4(y-7) && r.signature==signature;});
        if (rule==superdense_488_rules.end()) throw std::logic_error("Unknown superdense 4.8.8 boundary environment");
        if (rule->prepare_x) s.prepare_x.push_back(a);
        if (rule->measure_x) s.measure_x.push_back(a);
        for (const auto& gate:rule->gates)
            s.layers[gate[0]].emplace_back(at.at({x+gate[1],y+gate[2]}),at.at({x+gate[3],y+gate[4]}));
    }
    for (auto& layer:s.layers) std::sort(layer.begin(),layer.end());
    g.schedules.push_back(s);
    std::reverse(s.layers.begin(),s.layers.end());
    std::swap(s.prepare_x,s.measure_x);
    g.schedules.push_back(std::move(s));
    return g;
}

// Exact binary CSS Choi calculation. References are algebraic only; no extra
// physical qubits or feedback gates are introduced into the returned circuit.
struct Bits
{
    std::vector<uint64_t> words;
    explicit Bits(size_t n=0):words((n+63)/64) {}
    void flip(size_t i) {words.at(i/64)^=uint64_t{1}<<(i%64);}
    bool get(size_t i) const {return (words.at(i/64)>>(i%64))&1;}
    int top() const
    {
        for (size_t k=words.size(); k--;) if (words[k]) return 64*k+63-std::countl_zero(words[k]);
        return -1;
    }
    void add(const Bits& b) {for (size_t k=0;k<words.size();k++) words[k]^=b.words[k];}
    Qubits indices() const
    {
        Qubits result;
        for (size_t k=0;k<words.size()*64;k++) if (get(k)) result.push_back(k);
        return result;
    }
};
struct Row {Bits op,record;};
class Elimination
{
    std::vector<std::optional<Row>> pivots;
    size_t records;
public:
    std::vector<Bits> relations;
    Elimination(size_t n,size_t r):pivots(n),records(r) {}
    void insert(Bits v,Bits r)
    {
        while (v.top()>=0)
        {
            auto& p=pivots.at(v.top());
            if (!p) {p=Row{std::move(v),std::move(r)}; return;}
            v.add(p->op); r.add(p->record);
        }
        if (r.top()>=0) relations.push_back(std::move(r));
    }
    Bits query(Bits v) const
    {
        Bits r(records);
        while (v.top()>=0)
        {
            const auto& p=pivots.at(v.top());
            if (!p) throw std::logic_error("Superdense cycle does not extract a declared stabilizer/logical");
            v.add(p->op); r.add(p->record);
        }
        return r;
    }
    Bits reduce(Bits v) const
    {
        for (size_t b=pivots.size(); b--;) if (pivots[b] && v.get(b)) v.add(pivots[b]->op);
        return v;
    }
    std::vector<Bits> basis() const
    {
        std::vector<Bits> result;
        for (size_t b=pivots.size(); b--;) if (pivots[b]) result.push_back(pivots[b]->op);
        return result;
    }
};
struct Instrument
{
    std::array<std::vector<Bits>,2> input,output;
    std::array<Bits,2> logical;
    std::vector<std::pair<Bits,int>> flags;
};
Instrument instrument(const Geometry& g,const Schedule& s)
{
    const size_t n=g.data,N=g.coordinates.size(),na=N-n;
    Instrument maps;
    for (int axis=0;axis<2;axis++)
    {
        std::vector<Bits> rows;
        for (size_t q=0;q<n;q++)
        {
            Bits v(N+n); v.flip(q); v.flip(N+q); rows.push_back(std::move(v));
        }
        for (size_t q=n;q<N;q++)
            if ((std::find(s.prepare_x.begin(),s.prepare_x.end(),q)!=s.prepare_x.end())==(axis==0))
            {Bits v(N+n);v.flip(q);rows.push_back(std::move(v));}
        for (const auto& layer:s.layers) for (auto [a,b]:layer)
        {
            if (axis==1) std::swap(a,b);
            for (auto& v:rows) if (v.get(a)) v.flip(b);
        }
        Elimination e(N+n,na);
        for (auto& v:rows) e.insert(std::move(v),Bits(na));
        for (size_t q=n;q<N;q++)
            if ((std::find(s.measure_x.begin(),s.measure_x.end(),q)!=s.measure_x.end())==(axis==0))
            {Bits v(N+n),r(na);v.flip(q);r.flip(q-n);e.insert(std::move(v),std::move(r));}
        Elimination constraints(na,0);
        for (const auto& r:e.relations) constraints.insert(r,Bits());
        for (const auto& r:constraints.basis()) maps.flags.emplace_back(r,axis);
        for (const auto& f:g.faces)
        {
            Bits incoming(N+n),outgoing(N+n);
            for (auto q:f.support) {incoming.flip(N+q);outgoing.flip(q);}
            maps.input[axis].push_back(constraints.reduce(e.query(std::move(incoming))));
            maps.output[axis].push_back(constraints.reduce(e.query(std::move(outgoing))));
        }
        Bits logical(N+n);
        for (size_t q=0;q<n;q++) {logical.flip(q);logical.flip(N+q);}
        maps.logical[axis]=constraints.reduce(e.query(std::move(logical)));
    }
    return maps;
}

// Native 4.8.8 Bell-polarity optimization and data Hadamard frames, exactly
// as candidate._inputs in the source bundle; then recompute the instrument.
void native488(Geometry& g)
{
    auto& forward=g.schedules[0];
    std::set<int> feeds;
    for (auto [a,b]:forward.layers[1])
        if (std::find(g.routing.begin(),g.routing.end(),b)!=g.routing.end()) feeds.insert(a);
    std::set<uint32_t> prepare(forward.prepare_x.begin(),forward.prepare_x.end());
    for (auto& [a,b]:forward.layers[0]) if (feeds.contains(b))
    {
        prepare.erase(a); prepare.insert(b); std::swap(a,b);
    }
    forward.prepare_x.assign(prepare.begin(),prepare.end());
    auto reverse=forward;
    std::reverse(reverse.layers.begin(),reverse.layers.end());
    std::swap(reverse.prepare_x,reverse.measure_x);
    g.schedules[1]=std::move(reverse);
    for (auto& s:g.schedules)
    {
        std::vector<bool> frame(g.data);
        for (const auto& layer:s.layers) for (auto [a,b]:layer)
        {if (a<(int)g.data) frame[a]=false;if (b<(int)g.data) frame[b]=true;}
        for (size_t q=0;q<g.data;q++) if (frame[q]) s.output_h.push_back(q);
    }
    g.schedules[0].input_h=g.schedules[1].output_h;
    g.schedules[1].input_h=g.schedules[0].output_h;
}

struct Gate {std::string name; Qubits targets;};
using Moments=std::vector<std::vector<Gate>>;
Moments physical_cycle(const Geometry& g,const Schedule& s,bool native_cz)
{
    Qubits ancillas(g.coordinates.size()-g.data);
    std::iota(ancillas.begin(),ancillas.end(),g.data);
    Moments moments;
    auto append=[&](const std::string& name,const Qubits& targets)
    {if (!targets.empty()) moments.push_back({{name,targets}});};
    append("R",ancillas);
    if (!native_cz)
    {
        append("H",s.prepare_x);
        for (const auto& layer:s.layers)
        {
            Qubits pairs;
            for (auto [a,b]:layer) pairs.insert(pairs.end(),{uint32_t(a),uint32_t(b)});
            // Preserve even a physically all-idle entangling layer.
            moments.push_back(pairs.empty() ? std::vector<Gate>{}:std::vector<Gate>{{"CX",pairs}});
        }
        append("H",s.measure_x);
    }
    else
    {
        std::set<uint32_t> pending(s.prepare_x.begin(),s.prepare_x.end());
        auto toggle=[&](uint32_t q) {if (!pending.erase(q)) pending.insert(q);};
        for (auto q:s.input_h) toggle(q);
        for (const auto& layer:s.layers)
        {
            std::set<uint32_t> active,targets;
            Qubits rotations,pairs;
            for (auto [a,b]:layer) {active.insert(a);active.insert(b);targets.insert(b);pairs.insert(pairs.end(),{uint32_t(a),uint32_t(b)});}
            for (auto q:active) if (pending.contains(q)!=targets.contains(q)) rotations.push_back(q);
            append("H",rotations); append("CZ",pairs);
            for (auto q:active) pending.erase(q);
            pending.insert(targets.begin(),targets.end());
        }
        for (auto q:s.measure_x) toggle(q);
        for (auto q:s.output_h) toggle(q);
        append("H",Qubits(pending.begin(),pending.end()));
        // Schedule as early as each physical wire permits; measurements share
        // a final moment. This is candidate.compact, without a scheduling solve.
        Moments compact;
        std::vector<int> last(g.coordinates.size(),-1);
        for (const auto& moment:moments) for (const auto& gate:moment)
        {
            const size_t stride=gate.name=="CZ" ? 2:1;
            for (size_t j=0;j<gate.targets.size();j+=stride)
            {
                Qubits qs(gate.targets.begin()+j,gate.targets.begin()+j+stride);
                int depth=0;
                for (auto q:qs) depth=std::max(depth,last[q]+1);
                if (compact.size()<size_t(depth+1)) compact.resize(depth+1);
                compact[depth].push_back({gate.name,qs});
                for (auto q:qs) last[q]=depth;
            }
        }
        moments=std::move(compact);
    }
    append("M",ancillas);
    return moments;
}

void validate_geometry(const Geometry& g)
{
    if (g.faces.size()!=(g.data-1)/2) throw std::logic_error("Incorrect superdense stabilizer census");
    Elimination rank(g.data,0);
    for (size_t f=0;f<g.faces.size();f++)
    {
        Bits v(g.data);
        for (auto q:g.faces[f].support) v.flip(q);
        rank.insert(v,Bits());
        for (size_t k=0;k<f;k++)
        {
            size_t overlap=0;
            for (auto q:g.faces[k].support) overlap+=v.get(q);
            if (overlap%2 || (overlap && g.faces[f].color==g.faces[k].color))
                throw std::logic_error("Superdense faces do not commute or are not three-colored");
        }
    }
    if (rank.basis().size()!=g.faces.size()) throw std::logic_error("Dependent superdense stabilizers");
    for (const auto& s:g.schedules) for (const auto& layer:s.layers)
    {
        std::set<int> used;
        for (auto [a,b]:layer)
            if (a==b || !used.insert(a).second || !used.insert(b).second)
                throw std::logic_error("Superdense CNOT layer has conflicting operations");
    }
}

stim::Circuit memory(uint32_t distance,uint32_t rounds,const CircuitNoisePolicy& policy,
    bool memory_x,bool octagons,bool opposite,const ComponentErrorDistribution& noise,
    uint32_t start_phase,bool native_cz)
{
    if (distance<3 || distance%2==0 || !rounds || (octagons && distance%4!=3) ||
        start_phase>(octagons ? 1u:0u))
        throw std::invalid_argument("superdense: need rounds>=1, odd d>=3 (4.8.8: d=4m-1), valid start phase");
    policy.validate(rounds);
    if (policy.check_implementation()!=CheckImplementation::GateSequence)
        throw std::invalid_argument("superdense requires a gate-sequence noise policy; use ideal-check generators for capacity/pheno");
    auto g=octagons ? geometry488(distance):geometry666(distance);
    if (native_cz && octagons) native488(g);
    validate_geometry(g);
    const size_t n=g.data,N=g.coordinates.size(),na=N-n;
    Qubits data(n);std::iota(data.begin(),data.end(),0);
    std::vector<Instrument> maps;
    std::vector<Moments> cycles;
    for (const auto& s:g.schedules) {maps.push_back(instrument(g,s));cycles.push_back(physical_cycle(g,s,native_cz));}
    stim::Circuit c;
    for (size_t q=0;q<N;q++)
        c.safe_append_u("QUBIT_COORDS",{uint32_t(q)},{double(g.coordinates[q].first),double(g.coordinates[q].second),0,q<n ? 0.:1.});
    auto moment=[&](const std::vector<Gate>& gates)
    {
        std::set<uint32_t> busy,collapsed;
        for (const auto& gate:gates)
        {
            busy.insert(gate.targets.begin(),gate.targets.end());
            const auto kind=gate.name=="R" ? CircuitOperationKind::Reset : gate.name=="M" ? CircuitOperationKind::Measure :
                (gate.name=="CX" || gate.name=="CZ") ? CircuitOperationKind::Unitary2:CircuitOperationKind::Unitary1;
            policy.append(c,{kind,gate.targets,gate.name});
            if (gate.name=="R" || gate.name=="M") collapsed.insert(gate.targets.begin(),gate.targets.end());
            if (gate.name=="M") policy.append(c,{CircuitOperationKind::PostMeasurement,gate.targets});
        }
        Qubits idle;
        for (size_t q=0;q<N;q++) if (!busy.contains(q)) idle.push_back(q);
        if (!idle.empty())
        {
            policy.append(c,{CircuitOperationKind::GateIdle,idle});
        }
        // Packed native moments can reset one ancilla while other qubits
        // have gates. Every noncollapsed qubit gets the resonator-wait noise.
        if (!collapsed.empty())
        {
            Qubits waiting;
            for (size_t q=0;q<N;q++) if (!collapsed.contains(q)) waiting.push_back(q);
            if (!waiting.empty()) policy.append(c,{CircuitOperationKind::MeasurementIdle,waiting});
        }
        c.safe_append_u("TICK",{});
    };
    moment({{"R",data}});
    auto rotations=[&](const Qubits& frame)
    {
        std::set<uint32_t> h(frame.begin(),frame.end());
        if (memory_x) for (auto q:data) {if (!h.erase(q)) h.insert(q);}
        if (!h.empty()) moment({{"H",Qubits(h.begin(),h.end())}});
    };
    rotations(g.schedules[start_phase].input_h);
    uint64_t measured=0;
    auto append_records=[&](Qubits& records,const Bits& bits,uint64_t offset)
    {
        for (auto j:bits.indices()) records.push_back(stim::GateTarget::rec(-(int64_t(measured)-int64_t(offset+j))).data);
    };
    auto center=[&](const Face& f)
    {
        double x=0,y=0;
        for (auto q:f.support) {x+=double(g.coordinates[q].first)/f.support.size();y+=double(g.coordinates[q].second)/f.support.size();}
        return std::pair{x,y};
    };
    auto detector=[&](Qubits records,double x,double y,size_t rd,int color,int basis,bool flag)
    {c.safe_append_u("DETECTOR",records,{x,y,double(rd),double(color),double(basis),flag ? 1.:0.});};
    for (size_t rd=0;rd<rounds;rd++)
    {
        const size_t phase=(start_phase+rd)%maps.size();
        const auto& map=maps[phase];
        policy.append(c,{CircuitOperationKind::RoundStart,data});
        const auto offset=measured;
        for (const auto& m:cycles[phase]) moment(m);
        measured+=na;
        for (const auto& [flag,axis]:map.flags)
        {
            Qubits records;append_records(records,flag,offset);
            double x=0,y=0;
            auto support=flag.indices();
            for (auto j:support) {x+=double(g.coordinates[n+j].first)/support.size();y+=double(g.coordinates[n+j].second)/support.size();}
            // Flags are colorless. Label the hook-error basis, opposite to
            // the stabilizer sector detected by this ideal record constraint.
            detector(records,x,y,rd,-1,1-axis,true);
        }
        for (size_t f=0;f<g.faces.size();f++) for (int axis=0;axis<2;axis++)
        {
            if ((!rd || !opposite) && ((axis==0)!=memory_x)) continue;
            Qubits records;append_records(records,map.input[axis][f],offset);
            if (rd) append_records(records,maps[(phase+maps.size()-1)%maps.size()].output[axis][f],offset-na);
            const auto [x,y]=center(g.faces[f]);
            detector(records,x,y,rd,g.faces[f].color,axis,false);
        }
        Qubits records;append_records(records,map.logical[memory_x ? 0:1],offset);
        if (!records.empty()) c.safe_append_ua("OBSERVABLE_INCLUDE",records,0);
    }
    const size_t final_phase=(start_phase+rounds-1)%maps.size();
    rotations(g.schedules[final_phase].output_h);
    moment({{"M",data}});measured+=n;
    for (size_t f=0;f<g.faces.size();f++)
    {
        Qubits records;
        append_records(records,maps[final_phase].output[memory_x ? 0:1][f],measured-n-na);
        for (auto q:g.faces[f].support) records.push_back(stim::GateTarget::rec(-int64_t(n-q)).data);
        const auto [x,y]=center(g.faces[f]);
        detector(records,x,y,rounds,g.faces[f].color,memory_x ? 0:1,false);
    }
    Qubits logical;
    for (auto q:data) logical.push_back(stim::GateTarget::rec(-int64_t(n-q)).data);
    c.safe_append_ua("OBSERVABLE_INCLUDE",logical,0);
    return apply_component_error_distribution(c,noise);
}
}

stim::Circuit cc666_superdense(uint32_t d,uint32_t r,const CircuitNoisePolicy& policy,bool x,
    bool opposite,const ComponentErrorDistribution& noise,bool native_cz)
{return memory(d,r,policy,x,false,opposite,noise,0,native_cz);}
stim::Circuit cc488_superdense(uint32_t d,uint32_t r,const CircuitNoisePolicy& policy,bool x,
    bool opposite,const ComponentErrorDistribution& noise,uint32_t start_phase,bool native_cz)
{return memory(d,r,policy,x,true,opposite,noise,start_phase,native_cz);}
