using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ViewChage : MonoBehaviour
{
    // Start is called before the first frame update\
    public Camera Cam1, Cam2;
    void Start()
    {
        Cam1.enabled = true;
        Cam2.enabled = false;
    }

    // Update is called once per frame
    void Update()
    {
        if (Input.GetKeyDown(KeyCode.Q))
        {
            if(Cam1.enabled) { 
                Cam1.enabled = false;
                Cam2.enabled = true;
            }
            else
            {
                Cam1.enabled = true;
                Cam2.enabled = false;
            }
        }
    }
}
