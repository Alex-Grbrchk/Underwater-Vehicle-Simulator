using System.Collections;
using System.Collections.Generic;
using UnityEngine.UI;
using UnityEngine;
using TMPro;


public class GUI : MonoBehaviour
{
    public Button UpHeight, DownHeight, UpSpeed, DownSpeed;
    public Text CurH, CurS;
    public TMP_InputField H, S;
    public Moving Mvs;//moving script
    
    void Start()
    {
        UpHeight.onClick.AddListener(UPH);
        DownHeight.onClick.AddListener(DWH);
        UpSpeed.onClick.AddListener(UPS);
        DownSpeed.onClick.AddListener(DWS);
    }

    // Update is called once per frame
    void Update()
    {
        CurH.text = "" + Mvs.transform.position.y + " m";
        if(float.Parse(H.text) <= 0)
        {
            Mvs.height = float.Parse(H.text);
        }

        CurS.text = "" + System.Math.Sqrt(System.Math.Pow(Mvs.RB.velocity.x, 2) + System.Math.Pow(Mvs.RB.velocity.z, 2))  + " mps";
        if (float.Parse(S.text) <= 3 && float.Parse(S.text) >=0)
        {
            Mvs.speed = float.Parse(S.text);
        }
        
    }
    void UPH()
    {
        if (Mvs.height < 0) {
            Mvs.height++;
        }
        else {
            if(Mvs.height > 0) Mvs.height = 0; 
        }
        H.text = "" + Mvs.height;
    }

    void DWH()
    {
        Mvs.height--;
        H.text = "" + Mvs.height;
    }

    void UPS()
    {
        if (Mvs.speed < 3)
        {
            Mvs.speed = Mvs.speed + 0.1f;
        }
        S.text = "" + Mvs.speed;
    }

    void DWS()
    {
        if (Mvs.speed > 0)
        {
            Mvs.speed = Mvs.speed - 0.1f;
        }
        else
        {
            if (Mvs.speed < 0) Mvs.speed = 0f;
        }
        S.text = "" + Mvs.speed;

    }
}
