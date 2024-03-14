using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ViewChage : MonoBehaviour
{
    // Start is called before the first frame update\
    public GameObject Cam1, Cam2;
    public Button Front, Main;
    void Start()
    {
        Cam1.SetActive (false);
        Cam2.SetActive (true);
        Front.gameObject.SetActive(true);
        Main.gameObject.SetActive(false);
        Front.onClick.AddListener(ViewCh1);
        Main.onClick.AddListener(ViewCh2);
    }

    // Update is called once per frame
    void Update()
    {
        if (Input.GetKeyDown(KeyCode.Q))
        {
            if (Cam1.activeSelf == true)
            {
                Cam1.SetActive(false);
                Cam2.SetActive(true);
            }
            else
            {
                Cam1.SetActive(true);
                Cam2.SetActive(false);
            }
        }
    }

    void ViewCh1()
    {

        Cam1.SetActive(true);
        Cam2.SetActive(false);
        Main.gameObject.SetActive(true);
        Front.gameObject.SetActive(false);
    }
    void ViewCh2()
    {
        Cam1.SetActive(false);
        Cam2.SetActive(true);
        Front.gameObject.SetActive(true);
        Main.gameObject.SetActive(false);
    }
}
