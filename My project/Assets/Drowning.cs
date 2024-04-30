using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Numerics;
using UnityEngine;


/*      ===     Динамическая модель     ===
*/

public class DSRV
{
    public double reff; //глубина
    
    public string controlDescription;

    //string names;
    public double L;//длина
    public double deltaMax;//максимальный угол поворота
    public double T_delta;
    public double U0;//круизная скорость м/с
    public double W0;
    public List<double> nu;//вектор скорости
    public List<double> u_actual;

    //public Vector<string> controls;
    public int dimU;


    public double Iy;
    public double m;
    public double Mqdot;
    public double Zqdot;
    public double Mwdot;
    public double Zwdot;

    public double m11;
    public double m12;
    public double m22;
    public double m21;

    public double detM;


    public double Mq;
    public double Zq;
    public double Mw;
    public double Zw;
    public double Mdelta;
    public double Zdelta;


    public double e_int;
    public double wn;
    public double zeta;


    public double w_max;
    public double z_d;
    public double w_d;
    public double a_d;
    public double wn_d;
    public double zeta_d;





    public DSRV(double z_dN)
    {
        
        reff = z_dN;
       
        L = 5;
        deltaMax = 50;//20;
        T_delta =  Time.deltaTime;//1;
        U0 = 4.11;
        W0 = 0;

        nu = new List<double> { U0, 0, W0, 0, 0, 0 };
        u_actual = new List<double> { 0 };
        
        dimU = 1;

        Iy = 0.001925;
        m = 0.036391;
        Mqdot = -0.001573;
        Zqdot = -0.000130;
        Mwdot = -0.000146;
        Zwdot = -0.031545;

        m11 = m - Zwdot;
        m12 = -Zqdot;
        m22 = Iy - Mqdot;
        m21 = -Mwdot;

        detM = m11 * m22 - m12 * m21;


        Mq = -0.01131;
        Zq = -0.017455;
        Mw = 0.011175;
        Zw = -0.043938;
        Mdelta = -0.012797;
        Zdelta = 0.027695;


        e_int = 0;
        wn = 1;
        zeta = 1;


        w_max = 1;
        z_d = 0;
        w_d = 0;
        a_d = 0;
        wn_d = wn / 5;
        zeta_d = 1;
    }

    List<List<double>> Rzyx(in double phi, in double theta, in double psi)
    {
        double cphi = Math.Cos(phi);
        double sphi = Math.Sin(phi);
        double cth = Math.Cos(theta);
        double sth = Math.Sin(theta);
        double cpsi = Math.Cos(psi);
        double spsi = Math.Sin(psi);

        List<List<double>> R = new List<List<double>>
        {
            new() { cpsi * cth, -spsi * cphi + cpsi * sth * sphi, spsi * sphi + cpsi * cphi * sth },
            new() { spsi * cth, cpsi * cphi + sphi * sth * spsi, -cpsi * sphi + sth * spsi * cphi },
            new() { -sth, cth * sphi, cth * cphi }
        };

        return R;
    }
    List<List<double>> Tzyx(in double phi, in double theta)
    {
        double cphi = Math.Cos(phi);
        double sphi = Math.Sin(phi);
        double cth = Math.Cos(theta);
        double sth = Math.Sin(theta);
        List<List<double>> T = new List<List<double>>
        {
            new() { 1, sphi * sth / cth, cphi * sth / cth },
            new() { 0, cphi, -sphi },
            new() { 0, sphi / cth, cphi / cth }
        };
        return T;
    }

    List<List<double>> Matmul3(in List<List<double>> A, in List<List<double>> B)
    {
        List<List<double>> T2 = new();
        for (int i = 0; i < A.Count; i++)
        {
            List<double> T1 = new();
            for (int j = 0; j < B[0].Count; j++)
            {
                T1.Add(0);
                for (int k = 0; k < A[0].Count; k++)
                {
                    T1[j] = T1[j] + A[i][k] * B[k][j];
                }
            }
            T2.Add(T1);
        }

        return T2;
    }
    List<double> Matmul(in List<List<double>> A, in List<double> B)
    {

        List<List<double>> A1 = new();
        for (int i = 0; i < B.Count; i++)
        {
            List<double> T1 = new();
            for (int j = 0; j < A[0].Count; j++)
            {
                T1.Add(0);
            }
            A1.Add(T1);
        }
        for (int i = 0; i < B.Count; i++)
        {
            A1[i][0] = B[i];
        }
        List<List<double>> A2 = Matmul3(A, A1);

        List<double> MM = new();
        for (int i = 0; i < A2.Count; i++)
        {
            MM.Add(A2[i][0]);
        }
        return MM;
    }

    void AttitudeEuler(ref List<double> eta, ref List<double> nu, in double sampleTime)
    {
        List<double> U_nu = new(){ nu[0],nu[1],nu[2] };
        List<double> p_dot = Matmul(Rzyx(eta[3], eta[4], eta[5]), U_nu);
        List<double> U_nu2 = new() { nu[3],nu[4],nu[5] };
        List<double> v_dot = Matmul(Tzyx(eta[3], eta[4]), U_nu2);


        eta[0] = eta[0] + sampleTime * p_dot[0];
        eta[1] = eta[1] + sampleTime * p_dot[1];
        eta[2] = eta[2] + sampleTime * p_dot[2];
        eta[3] = eta[3] + sampleTime * v_dot[0];
        eta[4] = eta[4] + sampleTime * v_dot[1];
        eta[5] = eta[5] + sampleTime * v_dot[2];
    }

    public List<double> Simulate(in List<double> Data, in double N)
    {
                            
        double t = N;                      // Время на 1 такт

        //Векторы начального состояния
        List<double> eta = Data;    

        double u_control;
        u_control = DepthAutopilot(eta, nu, t);
            
        //передача динамики транспортного средства и ориентации
        Dynamics(eta,ref nu,ref u_actual, u_control, t);// пересчитывает S_u_actual и S_nu
        AttitudeEuler(ref eta,ref nu, t);//пересчитывает eta
       
        return eta;
    }

    void RefModel3(ref double  x_d, ref double v_d, ref double a_d, in double r, in double wn_d, in double zeta_d, in double v_max, in double sampleTime)
    {


        double j_d = Math.Pow(wn_d, 3) * (r - x_d) - (2 * zeta_d + 1) * Math.Pow(wn_d, 2) * v_d - (2 * zeta_d + 1) * wn_d * a_d;


        x_d = x_d + sampleTime * v_d;
        v_d = v_d + sampleTime * a_d;
        a_d = a_d + sampleTime * j_d;


        if (v_d > v_max)
        {
            v_d = v_max;
        }
        else if (v_d < -v_max)
        {
            v_d = -v_max;
        }
    }

    double PIDpolePlacement(
    ref double e_int,
    in double e_x,
    in double e_v,
    ref double x_d,
    ref double v_d,
    ref double a_d,
    in double m,
    in double d,
    in double k,
    in double wn_d,
    in double zeta_d,
    in double wn,
    in double zeta,
    in double r,
    in double v_max,
    in double sampleTime
    )

    {
        double Kp = m * Math.Pow(wn, 2) - k;
        double Kd = m * 2 * zeta * wn - d;
        double Ki = (wn / 10) * Kp;


        double u = -Kp * e_x - Kd * e_v - Ki * e_int;


        e_int = e_int + sampleTime * e_x;


        RefModel3(ref x_d, ref v_d,ref a_d, r, wn_d, zeta_d, v_max, sampleTime);

        return u;
    }
    double DepthAutopilot(in List<double> eta, in List<double> nu, in double sampleTime)
    {
        double delta_c = 0;

        double z = eta[2];
        double w = nu[2];
        double e_z = z - z_d;
        double e_w = w - w_d;
        double r = reff;


        double m = m11;
        double d = 0;
        double k = 0;


        delta_c = PIDpolePlacement(
            ref e_int,
            e_z,
            e_w,
            ref z_d,
            ref w_d,
            ref a_d,
            m,
            d,
            k,
            wn_d,
            zeta_d,
            wn,
            zeta,
            r,
            w_max,
            sampleTime
            );

        return delta_c;
    }
    void Dynamics(in List<double> eta, ref List<double> nu, ref List<double> u_actual, in double u_control, in double sampleTime)
    {
        double delta_c = u_control;
        double delta = u_actual[0];
        double w = nu[2];
        double q = nu[4];
        double theta = eta[4];


        double U = Math.Sqrt(Math.Pow(U0, 2) + Math.Pow((W0 + w), 2));


        double Mtheta = -0.156276 / Math.Pow(U, 2);


        if (Math.Abs(delta) >= deltaMax * Math.PI / 180)
        {
            delta = Math.Sign(delta) * deltaMax * Math.PI / 180;

        }


        double Z = Zq * q + Zw * w + Zdelta * delta;
        double M = Mq * q + Mw * w + Mtheta * theta + Mdelta * delta;


        List<double> nu_dot = new(){ 0,0,0,0,0,0 };
        nu_dot[2] = (m22 * Z - m12 * M) / detM;
        nu_dot[4] = (-m21 * Z + m11 * M) / detM;


        double delta_dot = (delta_c - delta) / T_delta;


        nu[0] = nu[0] + sampleTime * nu_dot[0];
        nu[1] = nu[1] + sampleTime * nu_dot[1];
        nu[2] = nu[2] + sampleTime * nu_dot[2];
        nu[3] = nu[3] + sampleTime * nu_dot[3];
        nu[4] = nu[4] + sampleTime * nu_dot[4];
        nu[5] = nu[5] + sampleTime * nu_dot[5];
        delta = delta + sampleTime * delta_dot;


        nu[0] = U0;

        u_actual[0] = delta;
    }
};


public class Drowning : MonoBehaviour
{
    //public Rigidbody RB;
    public GameObject Boat;
    DSRV Vehicle;
    List<Double> BoatD;
    List<Double> BoatG;
    public double height;

    void Start()
    {
        //RB = gameObject.GetComponent<Rigidbody>();
        height = 0;
    }
    private void Awake()
    {
        
        Vehicle = new(Boat.transform.position.y);
        BoatD = new() { Boat.transform.position.x, Boat.transform.position.z, Boat.transform.position.y, Boat.transform.rotation.x*Math.PI/180, Boat.transform.rotation.z * Math.PI / 180, Boat.transform.rotation.y * Math.PI / 180 };
        BoatG = BoatD;
        Vehicle.U0 = 0.0001; //скорость
    }
    bool Sim = false;
    private void FixedUpdate()
    {
        if (Input.GetKey(KeyCode.W)) {

            if (Boat.transform.position.y != height)
            {
                Vehicle.reff = height;
            }
            Vehicle.U0 = 4.11;
            Sim = true;
        }
        if (Sim)
        {
            BoatG = Vehicle.Simulate(BoatD, Time.fixedDeltaTime);
            BoatD = BoatG;
            Boat.GetComponent<Rigidbody>().velocity = new UnityEngine.Vector3((float)BoatG[1], (float)BoatG[2], (float)BoatG[0]) - gameObject.transform.position;
            Boat.GetComponent<Rigidbody>().AddTorque(UnityEngine.Quaternion.ToEulerAngles(UnityEngine.Quaternion.Euler((float)(BoatG[4] * 180 / Math.PI), (float)(BoatG[5] * 180 / Math.PI), (float)(BoatG[3] * 180 / Math.PI))) - UnityEngine.Quaternion.ToEulerAngles(gameObject.transform.rotation), ForceMode.VelocityChange);
        }

        //Boat.GetComponent<Rigidbody>().Move(new UnityEngine.Vector3((float)BoatG[1], (float)BoatG[2], (float)BoatG[0]), UnityEngine.Quaternion.Euler((float)(-BoatG[3] * 180 / Math.PI), (float)(-BoatG[5] * 180 / Math.PI), (float)(-BoatG[4] * 180 / Math.PI)));
        //Boat.transform.SetPositionAndRotation(new UnityEngine.Vector3((float)BoatG[1], (float)BoatG[2], (float)BoatG[0]), UnityEngine.Quaternion.Euler((float)(-BoatG[3]*180/Math.PI ), (float)( -BoatG[5] * 180 / Math.PI), (float)( -BoatG[4] * 180 / Math.PI)));
        //BoatD = new() { Boat.transform.position.x, Boat.transform.position.z, Boat.transform.position.y, Boat.transform.rotation.x * Math.PI / 180, Boat.transform.rotation.z * Math.PI / 180, Boat.transform.rotation.y * Math.PI / 180 };
    }
    private void OnCollisionEnter(Collision collision)
    {
        Sim = false;
        BoatD = new() { Boat.transform.position.x, Boat.transform.position.z, Boat.transform.position.y, Boat.transform.rotation.x * Math.PI / 180, Boat.transform.rotation.z * Math.PI / 180, Boat.transform.rotation.y * Math.PI / 180 };
    }
}






/*     ===      Кинематическая модель       ===
*/



    /*double[] xt = new double[6];
    double[] zadt = new double[17];
    double t;

    public const double pi = 3.141592653589793;           // число pi
    public const double pi2 = 6.283185307179586;          // число 2*pi
    public const double pi05 = 1.5707963267948966;        // число pi/2
    public const double rad_gr = 57.295779513;            // перевод радиан в градусы
    public const double gr_rad = 0.017453292519943;       // перевод градусов в радианы
    public const double yz_m = 0.514;                 // перевод узлов в м/c
    public const double m_yz = 1.9455252;                // перевод м/с в узлы

    public double d_t = 2.25;

    // контроль принадлежности угла диапазону [0,2*pi]
    double c_ygol(double a)
    {
        if (a < 0) return a += pi2;
        if (a > pi2) return a -= pi2;
        if (Math.Abs(a) < 1e-6) return 0;
        return a;
    }

    // определение заданного дифферента при маневрe по глубине
    void different(ref double[] x, ref double[] z)
    {
        double R1 = z[5], R2 = z[6], tetp = z[8], d_h = Math.Abs(x[2] - z[0]);
        if (d_h > R1 + R2) { z[3] = tetp; return; }
        double tet = Math.Acos(1 - d_h / (R1 + R2));
        z[3] = (tet > tetp) ? tetp : tet;
    }

    // расчет кинематических параметров
    // x - вектор кинематических параметров
    // z - вектор заданных параметров движения
    // zme - учитывается возможность движения по змейке (0 - нет, 1 - да)
    void kinem(ref double[] x, ref double[] z, int zme)
    { // маневр по скорости
        if (Math.Abs(z[1] - x[3]) > 1e-3)
        {
            double tv = Math.Abs(z[1] - x[3]) / z[7];
            double aV = (z[1] > x[3]) ? z[7] : -z[7];
            double A, V = x[3], zV = z[1];
            if (tv < d_t)
            { A = V * tv + 0.5f * aV * tv * tv + zV * (d_t - tv); x[3] = zV; }
            else
            { A = V * d_t + 0.5f * aV * d_t * d_t; x[3] += aV * d_t; }
            x[0] += A * Math.Cos(x[4]); x[1] += A * Math.Sin(x[4]);
            return;
        }
        // маневр по курсу
        if (Math.Abs(z[2] - x[4]) > 1e-3)
        {
            double wK, f, fs, a, e, tk;
            double swK = x[3] / z[4];
            if (z[2] > x[4])
            {
                if (z[2] - x[4] <= pi)
                { wK = swK; tk = Math.Abs((z[2] - x[4]) / wK); }
                else
                { wK = -swK; tk = Math.Abs((x[4] - z[2] + pi2) / wK); }
            }
            else
            {
                if (x[4] - z[2] <= pi)
                { wK = -swK; tk = Math.Abs((z[2] - x[4]) / wK); }
                else
                { wK = swK; tk = Math.Abs((z[2] - x[4] + pi2) / wK); }
            }
            if (tk > d_t)
            {
                f = swK * d_t; fs = 0.5 * f * wK / swK; a = 2 * z[4] * Math.Sin(0.5 * f);
                x[0] += a * Math.Cos(x[4] + fs); x[1] += a * Math.Sin(x[4] + fs);
                x[4] = c_ygol(x[4] += 2 * fs);
            }
            else
            {
                f = swK * tk; fs = 0.5f * f * wK / swK; a = 2 * z[4] * Math.Sin(0.5 * f);
                e = x[4] + fs;
                x[0] += a * Math.Cos(e) + x[3] * Math.Cos(z[2]) * (d_t - tk);
                x[1] += a * Math.Sin(e) + x[3] * Math.Sin(z[2]) * (d_t - tk);
                x[4] = c_ygol(x[4] += 2 * fs);
            }
            return;
        }
        // маневр по глубине
        if (Math.Abs(z[0] - x[2]) > 1e-2)
        {
            double R1 = z[5], R2 = z[6], hz = z[0], h = x[2], tet = x[5], tetz = z[3];
            double wt, tt, tp, tz, A, B;
            if (h < hz)  // погружение
            {
                if (tet < tetz - 0.001)
                {
                    if (h < hz - R2 * (1 - Math.Cos(tetz)))
                    {
                        wt = x[3] / R1; tt = Math.Abs((tetz - tet) / wt);
                        if (tt > d_t)
                        {
                            A = R1 * (Math.Sin(tet + wt * d_t) - Math.Sin(tet));
                            B = R1 * (Math.Cos(tet) - Math.Cos(tet + wt * d_t));
                            x[5] += wt * d_t;
                        }
                        else
                        {
                            B = R1 * (Math.Cos(tet) - Math.Cos(tet + wt * tt));
                            if (h + B + x[3] * (d_t - tt) * Math.Sin(tetz) > hz - R2 * (1 - Math.Cos(tetz)))
                            {
                                if (Math.Abs(hz - h - B - R2 * (1 - Math.Cos(tetz))) < 0.01)
                                {
                                    wt = x[3] / R1; tz = tetz / wt;
                                    A = R1 * (Math.Sin(tet + wt * tt) - Math.Sin(tet)) +
                                        x[3] * (d_t - tt - tz) + R2 * (Math.Sin(tet + wt * tz) - Math.Sin(tet));
                                    B += R2 * (Math.Cos(tet) - Math.Cos(tet + wt * tz));
                                    x[5] = tet;
                                }
                                else
                                {
                                    wt = x[3] / R2; tz = tetz / wt;
                                    tp = (hz - h - B - R2 * (1 - Math.Cos(tetz))) / (x[3] * Math.Sin(tetz));
                                    if (d_t - tt - tp < tz)
                                    {
                                        A = R1 * (Math.Sin(tet + wt * tt) - Math.Sin(tet)) +
                                          x[3] * tp * Math.Cos(tetz) +
                                          R2 * (Math.Sin(tetz) - Math.Sin(wt * (d_t - tt - tp)));
                                        B += x[3] * tp * Math.Sin(tetz) +
                                                     R2 * (Math.Cos(tetz) - Math.Cos(wt * (d_t - tt - tp)));
                                        x[5] = tetz - wt * (d_t - tt - tp);
                                    }
                                    else
                                    {
                                        A = R1 * (Math.Sin(tet + wt * tt) - Math.Sin(tet)) +
                                          x[3] * tp * Math.Cos(tetz) + x[3] * (d_t - tt - tp - tz) +
                                          R2 * (Math.Sin(tetz) - Math.Sin(wt * tz));
                                        B += x[3] * tp * Math.Sin(tetz) + R2 * (1 - Math.Cos(tetz));
                                        x[5] = 0;
                                    }
                                }
                            }
                            else
                            {
                                A = R1 * (Math.Sin(tet + wt * tt) - Math.Sin(tet)) + x[3] * (d_t - tt) * Math.Cos(tetz);
                                B = R1 * (Math.Cos(tet) - Math.Cos(tet + wt * tt)) + x[3] * (d_t - tt) * Math.Sin(tetz);
                                x[5] += wt * tt;
                            }
                        }
                        x[0] += A * Math.Cos(x[4]); x[1] += A * Math.Sin(x[4]); x[2] += B;
                        return;
                    }
                }

                if (Math.Abs(tet - tetz) < 0.001)
                {
                    if (h + x[3] * Math.Sin(tetz) * d_t < hz - R2 * (1 - Math.Cos(tetz)))
                    {
                        x[0] += x[3] * Math.Cos(tetz) * Math.Cos(x[4]) * d_t;
                        x[1] += x[3] * Math.Cos(tetz) * Math.Sin(x[4]) * d_t;
                        x[2] += x[3] * Math.Sin(tetz) * d_t;
                        return;
                    }
                    else
                    {
                        wt = x[3] / R2; tt = Math.Abs(tet / wt);
                        tz = Math.Abs((hz - h - R2 * (1 - Math.Cos(tetz))) / (x[3] * Math.Sin(tetz)));
                        if (tz + tt < d_t)
                        {
                            A = x[3] * Math.Cos(tetz) * tz + R2 * Math.Sin(tetz) + x[3] * (d_t - tt - tz);
                            B = x[3] * Math.Sin(tetz) * tz + R2 * (1 - Math.Cos(tetz));
                            x[0] += A * Math.Cos(x[4]); x[1] += A * Math.Sin(x[4]);
                            x[2] += B;
                            x[5] -= wt * tt;
                            return;
                        }
                        else
                        {
                            A = x[3] * Math.Cos(tetz) * tz + R2 * Math.Sin(tet) - R2 * Math.Sin(tet - wt * (d_t - tz));
                            B = x[3] * Math.Sin(tetz) * tz + R2 * Math.Cos(tet) - R2 * Math.Cos(tet - wt * (d_t - tz));
                            x[0] += A * Math.Cos(x[4]); x[1] += A * Math.Sin(x[4]);
                            x[2] += B;
                            x[5] -= wt * (d_t - tz);
                            return;
                        }
                    }
                }
            }
            else  // всплытие
            {
                if (tet < tetz - 0.001)
                {
                    if (h > hz + R1 * (1 - Math.Cos(tetz)))
                    {
                        wt = x[3] / R2; tt = Math.Abs((tetz - tet) / wt);
                        if (tt > d_t)
                        {
                            A = R2 * Math.Sin(tet + wt * d_t) - R2 * Math.Sin(tet);
                            B = R2 * Math.Cos(tet) - R2 * Math.Cos(tet + wt * d_t);
                            x[5] += wt * d_t;
                        }
                        else
                        {
                            B = R2 * (Math.Cos(tet) - Math.Cos(tet + wt * tt));
                            if (h - B - x[3] * (d_t - tt) * Math.Sin(tetz) < hz + R1 * (1 - Math.Cos(tetz)))
                            {
                                if (Math.Abs(hz - h + B + R1 * (1 - Math.Cos(tetz))) < 0.01)
                                {
                                    double wt2 = x[3] / R2, wt1 = x[3] / R1;
                                    tt = Math.Abs(tetz - tet) / wt2; tz = tetz / wt;
                                    A = R2 * (Math.Sin(tet + wt2 * tt) - Math.Sin(tet)) +
                                        x[3] * (d_t - tt - tz) + R1 * (Math.Sin(tet + wt * tz) - Math.Sin(tet));
                                    B += R1 * (Math.Cos(tet) - Math.Cos(tet + wt * tz));
                                    x[5] = 0;
                                }
                                else
                                {
                                    wt = x[3] / R1; tz = tetz / wt;
                                    tp = ((h - B - R1 * (1 - Math.Cos(tetz))) - hz) / (x[3] * Math.Sin(tetz));
                                    if (d_t - tt - tp < tz)
                                    {
                                        A = R2 * (Math.Sin(tet + wt * tt) - Math.Sin(tet)) +
                                          x[3] * tp * Math.Cos(tetz) +
                                          R1 * (Math.Sin(tetz) - Math.Sin(wt * (d_t - tt - tp)));
                                        B += x[3] * tp * Math.Sin(tetz) +
                                                     R1 * (Math.Cos(tetz) - Math.Cos(wt * (d_t - tt - tp)));
                                        x[5] = wt * (d_t - tt - tp);
                                    }
                                    else
                                    {
                                        A = R2 * (Math.Sin(tet + wt * tt) - Math.Sin(tet)) +
                                          x[3] * tp * Math.Cos(tetz) +
                                          x[3] * (d_t - tt - tp - tz) +
                                          R1 * (Math.Sin(tetz) - Math.Sin(wt * tz));
                                        B += x[3] * tp * Math.Sin(tetz) + R1 * (1 - Math.Cos(tetz));
                                        x[5] = 0;
                                    }
                                }
                            }
                            else
                            {
                                A = R2 * (Math.Sin(tet + wt * tt) - Math.Sin(tet)) +
                                  x[3] * (d_t - tt) * Math.Cos(tetz);
                                B = R2 * (Math.Cos(tet) - Math.Cos(tet + wt * tt)) +
                                    x[3] * (d_t - tt) * Math.Sin(tetz);
                                x[5] += wt * tt;
                            }
                        }
                        x[0] += A * Math.Cos(x[4]); x[1] += A * Math.Sin(x[4]); x[2] -= B;
                        return;
                    }
                    else
                    {
                        wt = x[3] / R1; tt = Math.Abs((tetz - tet) / wt);
                        if (tt > d_t)
                        {
                            A = R1 * Math.Sin(tet) - R1 * Math.Sin(tet - wt * d_t);
                            B = R1 * Math.Cos(tet - wt * d_t) - R1 * Math.Cos(tet);
                            x[5] -= wt * d_t;
                        }
                        else
                        {
                            A = R1 * Math.Sin(tet) - R1 * Math.Sin(tet - wt * tt) + x[3] * (d_t - tt);
                            B = R1 * Math.Cos(tet - wt * tt) - R1 * Math.Cos(tet);
                            x[5] = 0;
                        }
                        x[0] += A * Math.Cos(x[4]); x[1] += A * Math.Sin(x[4]); x[2] -= B;
                        return;
                    }
                }
                if (Math.Abs(tet - tetz) < 0.001)
                {
                    if (h - x[3] * Math.Sin(tetz) * d_t > hz + R1 * (1 - Math.Cos(tetz)))
                    {
                        x[0] += x[3] * Math.Cos(tetz) * Math.Cos(x[4]) * d_t;
                        x[1] += x[3] * Math.Cos(tetz) * Math.Sin(x[4]) * d_t;
                        x[2] -= x[3] * Math.Sin(tetz) * d_t;
                        return;
                    }
                    else
                    {
                        wt = x[3] / R1; tt = Math.Abs(tet / wt);
                        tz = Math.Abs((h - hz - R1 * (1 - Math.Cos(tetz))) / (x[3] * Math.Sin(tetz)));
                        if (tz + tt < d_t)
                        {
                            A = x[3] * Math.Cos(tetz) * tz + R1 * Math.Sin(tetz) + x[3] * (d_t - tt - tz);
                            B = x[3] * Math.Sin(tetz) * tz + R1 * (1 - Math.Cos(tetz));
                            x[0] += A * Math.Cos(x[4]); x[1] += A * Math.Sin(x[4]);
                            x[2] -= B;
                            x[5] -= wt * tt;
                            return;
                        }
                        else
                        {
                            A = x[3] * Math.Cos(tetz) * tz + R1 * Math.Sin(tet) - R1 * Math.Sin(tet - wt * (d_t - tz));
                            B = x[3] * Math.Sin(tetz) * tz + R1 * Math.Cos(tet) - R1 * Math.Cos(tet - wt * (d_t - tz));
                            x[0] += A * Math.Cos(x[4]); x[1] += A * Math.Sin(x[4]);
                            x[2] -= B;
                            x[5] -= wt * (d_t - tz);
                            return;
                        }
                    }
                }
            }
        }
        else
        { x[2] = z[0]; z[3] = 0; }
        if (zme == 1)
        { // движение по "змейке"
            double Kv = z[4] * Math.Sin(z[9] * z[10]) / (z[1] * z[10]);
            x[0] += Kv * z[1] * Math.Cos(x[4]) * d_t;
            x[1] += Kv * z[1] * Math.Sin(x[4]) * d_t;
            return;
        }
        // равномерное прямолинейное движение
        x[0] += x[3] * Math.Cos(x[4]) * d_t; x[1] += x[3] * Math.Sin(x[4]) * d_t;
    }
    
    // Start is called before the first frame update

    public GameObject Boat;
    void Awake()
    {
        // телеуправляемый подводный объект - ТПО
        // носитель
        xt[0] = Boat.transform.position.x;//0;              //* X (Nord)   ** не вводится **
        xt[1] = Boat.transform.position.z;//0;              //* Y          ** не вводится **
        xt[2] = Boat.transform.position.y;//0;              //* Z (вниз)
        xt[3] = 40 * yz_m;         //* V
        xt[4] = 10 * gr_rad;      //* K
        xt[5] = 0;              //* Дифферент  ** не вводится **
        zadt[0] = xt[0];    //* Глубина             ** не вводится **
        zadt[1] = xt[3];    //* Скорость            ** не вводится **
        zadt[2] = xt[4];    //* Курс                ** не вводится **
        zadt[3] = 0;         //* Дифферент           ** не вводится **
        zadt[4] = 21;        //* радиус циркуляции в горизонтальной плоскости
        zadt[5] = 21;       //* радиус циркуляции в вертикальной плоскости (малый)
        zadt[6] = 30;       //* радиус циркуляции в вертикальной плоскости (большой)
        zadt[7] = 7;        //* ускорение разгона торможения
        zadt[8] = 30 * gr_rad;//* предельный дифферент
        zadt[9] = 30 * gr_rad;//* угловая скрость циркуляции на "змейке"
        zadt[10] = 1.3333f;  //* полупериод "змейки"
        zadt[11] = 1000;      //* радиус реагирования ССН
        zadt[12] = 20 * gr_rad; //* половина угла раствора (горизонтальная плоскость)
        zadt[13] = 40 * gr_rad; //* половина угла раствора (вертикальная плоскость)
        zadt[14] = 40 * yz_m;   //* скорость V1
        zadt[15] = 24 * yz_m;   //* скорость V2
        zadt[16] = 600;       //* минимальный радиус реагирования ССН
    }

    bool MV = false;
    }
    private void Update()
    {
        if (Input.GetKeyDown(KeyCode.W))
        {
            zadt[0] = -50;
            different(ref xt, ref zadt);
            MV = true;
        }
    }
    private void FixedUpdate()
    {
        if (MV)
        {
            Boat.transform.position = new Vector3((float)xt[0], (float)xt[2], (float)xt[1]);
            Boat.transform.rotation = Quaternion.Euler(0, (float)(xt[4] * rad_gr), (float)(-xt[5] * rad_gr));
            d_t = Time.fixedDeltaTime;
            kinem(ref xt, ref zadt, 0);
        }
        
    } */