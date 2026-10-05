#include "core/navigation/ins_navigation.h"
#include "core/navigation/ins_cloud.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

using namespace lightning;
using namespace lightning::ins;
void Require(bool ok,const char* why) { if(!ok) { std::cerr<<"FAIL: "<<why<<'\n'; std::exit(1); } }
int main() {
    GeoReference equator(0,0,0);
    Require(equator.Forward(0,0,0).norm()<1e-9,"origin is zero");
    Require((equator.Forward(0,0,10)-Vec3d(0,0,10)).norm()<1e-8,"ellipsoid normal is up");
    const double angle=1e-5, r=6378137;
    const Vec3d analytic(r*std::sin(angle),0,r*(std::cos(angle)-1));
    Require((equator.Forward(0,angle*180/M_PI,0)-analytic).norm()<1e-7,"independent equatorial ECEF case");
    Mat3d far_rotation;
    equator.Forward(0,90,0,&far_rotation);
    Require((far_rotation*Vec3d::UnitX()+Vec3d::UnitZ()).norm()<1e-12,"local east rotates into negative fixed up at 90deg longitude");
    Require((far_rotation*Vec3d::UnitZ()-Vec3d::UnitX()).norm()<1e-12,"local up rotates into fixed east, not transpose");
    GeoReference site(31.2,121.4,20);
    Mat3d rotation;
    const Vec3d p=site.Forward(31.201,121.402,26,&rotation), llh=site.Reverse(p);
    Require((llh-Vec3d(31.201,121.402,26)).norm()<1e-7,"LLH roundtrip");
    Require((rotation.transpose()*rotation-Mat3d::Identity()).norm()<1e-12,"ENU vector rotation orthogonal");
    Require((VehicleAttitude(0,0,0)*Vec3d::UnitX()-Vec3d::UnitY()).norm()<1e-12,"north heading");
    Require((VehicleAttitude(90,0,0)*Vec3d::UnitX()-Vec3d::UnitX()).norm()<1e-12,"east heading");
    Require((VehicleAttitude(180,0,0)*Vec3d::UnitX()+Vec3d::UnitY()).norm()<1e-12,"south heading");
    Require((VehicleAttitude(270,0,0)*Vec3d::UnitX()+Vec3d::UnitX()).norm()<1e-12,"west heading");
    Require((VehicleAttitude(90,20,0)*Vec3d::UnitX()).z()>0,"nose up is up");
    Require((VehicleAttitude(90,0,20)*Vec3d::UnitY()).z()>0,"right down means left up");
    Sample s; s.time.fill(10);s.arrival.fill(20);s.valid.fill(true);
    s.system_state=2;s.satellite_status=4;s.heading=0;s.velocity=Vec3d(0,-2,0);
    QualityPolicy policy;
    Require(CheckSample(s,policy,10.01,20.01).empty(),"good sample");
    Require(std::abs(Convert(s,equator).SignedSpeed()+2)<1e-12,"signed reverse speed");
    for(int status:{0,1,2,3,5,6,7,8,9}) { auto a=s;a.satellite_status=status;Require(!CheckSample(a,policy,10.01,20.01).empty(),"reject nonfixed or heading invalid"); }
    for(int state:{0,1,3}) { auto a=s;a.system_state=state;Require(!CheckSample(a,policy,10.01,20.01).empty(),"reject reference-changing mode"); }
    for(std::size_t i=0;i<FieldCount;++i) { auto a=s;a.valid[i]=false;Require(!CheckSample(a,policy,10.01,20.01).empty(),"every field required"); }
    Require(!CheckSample(s,policy,10.01,21).empty(),"wall timeout while ROS clock paused");
    Require(!CheckSample(s,policy,11,20.01).empty(),"old receive stamps");
    Require(!CheckSample(s,policy,9,20.01).empty(),"future receive stamps");
    auto bad=s;bad.time[Longitude]-=0.03;Require(CheckSample(bad,policy,10.01,20.01)=="field_time_skew","cross epoch skew");
    bad=s;bad.position_sigma.x()=std::numeric_limits<double>::quiet_NaN();Require(CheckSample(bad,policy,10.01,20.01)=="invalid_sigma","NaN sigma");
    bad=s;bad.position_sigma.y()=0.101;Require(CheckSample(bad,policy,10.01,20.01)=="position_sigma","horizontal threshold");
    bad=s;bad.attitude_sigma.x()=0.51;Require(CheckSample(bad,policy,10.01,20.01)=="attitude_sigma","heading threshold");
    RecoveryGate gate(3);
    Require(!gate.Observe(1,true)&&!gate.Observe(2,true)&&gate.Observe(3,true),"recovery counts new samples");
    gate.Reject();Require(!gate.Good()&&!gate.Observe(4,true),"loss immediate recovery reset");
    Require(!gate.Observe(4,true),"duplicate not recovery progress");
    PoseBuffer buf(2,0.11);
    Pose a; a.stamp=1;a.global_valid=true;buf.Add(a);
    Pose b=a;b.stamp=1.1;b.map_from_rear=SE3(SO3(),Vec3d(1,0,0));buf.Add(b);
    Pose middle;Require(buf.At(1.05,middle)&&std::abs(middle.map_from_rear.translation().x()-0.5)<1e-12,"interpolation");
    Require(!buf.At(1.11,middle)&&!buf.At(0.99,middle),"no extrapolation");
    CloudPtr points(new PointCloudType);PointType v;v.x=9;v.y=0;v.z=0;v.time=0;points->push_back(v);
    v.x=8;v.time=100;points->push_back(v);CloudPtr out;SE3 end;
    Require(Deskew(points,1,1.1,SE3(SO3(),Vec3d(1,0,0)),buf,out,end),"deskew with sensor lever arm");
    Require(std::abs(out->points[0].x-9)<1e-6&&std::abs(out->points[1].x-9)<1e-6,"static wall restored in end rear frame");
    PoseBuffer turn(2,0.11);
    Pose start;start.stamp=2;start.global_valid=true;turn.Add(start);
    Pose finish=start;finish.stamp=2.1;
    finish.map_from_rear=SE3(SO3(Quatd(Eigen::AngleAxisd(M_PI/2,Vec3d::UnitZ()))),Vec3d::Zero());
    turn.Add(finish);
    CloudPtr rotating(new PointCloudType);
    v.x=9;v.y=0;v.time=0;rotating->push_back(v);
    v.x=-1;v.y=-10;v.time=100;rotating->push_back(v);
    Require(Deskew(rotating,2,2.1,SE3(SO3(),Vec3d(1,0,0)),turn,out,end),"rotating sensor lever arm deskew");
    for(const auto& point:*out) Require(std::abs(point.x)<1e-5&&std::abs(point.y+10)<1e-5,"ninety-degree turn restores stationary landmark");
    Require(!Deskew(points,1,1.2,SE3(),buf,out,end),"uncovered scan rejected");
    b.stamp=1.4;buf.Add(b);Require(!buf.Covers(1,1.4),"interior gap rejected");
    buf.Clear();Require(!buf.At(1,middle),"quality reset clears motion history");
    PoseBuffer wrap(2,1);a.stamp=1;a.map_from_rear=SE3(VehicleAttitude(359,0,0),Vec3d::Zero());
    b=a;b.stamp=2;b.map_from_rear=SE3(VehicleAttitude(1,0,0),Vec3d::Zero());wrap.Add(a);wrap.Add(b);
    Require(wrap.At(1.5,middle)&&(middle.map_from_rear.so3()*Vec3d::UnitX()-Vec3d::UnitY()).norm()<1e-10,"short quaternion path across north");
    std::cout<<"ins_navigation_test passed\n";
}
